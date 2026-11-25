#!/lab/bin/perl
#
#   snapshot.pl - making a snapshot of the PDB with today's date.
#
#   Copyright (C) 2001-2026 Degui Zhi and John-Marc Chandonia
#
#   This program is free software; you can redistribute it and/or modify it
#   under the terms of the GNU General Public License as published by the
#   Free Software Foundation; either version 2, or (at your option) any
#   later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
#
#
# This program is part of the PAST (PDB Archival Snapshot Toolkit) package.
#
# Contact:
# scope@compbio.berkeley.edu
#
# See README.txt for the options, the layout of a snapshot and the
# format of the manifest.
#
# Rule kept throughout: content that is already in the archive is never
# stored again and never gets a new modification time (ParsePDB makes a
# new release from either), and anything that would fetch again, or list
# as removed, more files than the limits allow stops the run before a
# snapshot is published.

use strict;
use warnings;
use LWP::UserAgent;
use HTTP::Request;
use HTTP::Date qw(str2time);
use Digest::MD5;
use Digest::SHA;
use Errno;
use File::Path qw(make_path remove_tree);
use File::Temp qw(tempdir);
use File::Copy qw(copy);
use File::Basename qw(dirname basename);
use Cwd qw(realpath);
use Getopt::Long;
use Fcntl qw(:flock);
use IO::Handle;
use POSIX qw(strftime);
use Time::HiRes ();
use Time::Local qw(timegm);

# ################## defaults ###########################
# Each can be set by the environment variable named here, or by the
# command-line option (see usage).

my %opt = (
    root             => $ENV{PAST_ROOT}       // "/lab/db/pdb",
    source           => $ENV{PAST_SOURCE}     // "https://files.wwpdb.org/pub/pdb",
    tmpdir           => $ENV{PAST_TMPDIR}     // "/tmp",
    onhold_url       => $ENV{PAST_ONHOLD_URL} // "https://www.rcsb.org/search/fastaDownload?pdbid=&authorList=&title=&status=&sequenceAvailable=&date=&excludeWDRN=yes&displayType=fasta",
    seqres_path      => "derived_data/pdb_seqres.txt.gz",
    index_path       => "derived_data/index/",
    index_required   => "author.idx,cmpd_res.idx,compound.idx,crystal.idx,entries.idx,on_hold.list,pending_waiting.list,resolu.idx,source.idx,status_query.csv,status_query.seq",
    gzip             => "gzip",
    date             => $ENV{PAST_DATE},
    # guards checked before anything is created (per kind of file)
    max_count_change => 0.10,   # change in the number of files, fraction of the previous number
    max_deletions    => 60,     # files that disappear
    max_candidates   => 0.20,   # files to check with HEAD, fraction of the listed number
    min_pdb_share    => 0.90,   # current entries that have a PDB-format file
    max_still_listed => 0,      # removed files whose entry is still in the same holdings file
    max_vanished     => 0,      # removed files whose entry is in neither holdings file
    max_reappeared   => 0,      # new files that an earlier snapshot listed as removed
    # limits checked during the run
    max_changed      => 0.10,   # files stored with changed content, fraction of the previous number
    max_repaired     => 5,      # files missing from the archive and fetched again (total)
    max_failures     => 20,     # failed files before the run stops
    max_compared     => 100,    # files downloaded and found to be the archived copy
    max_hours        => 144,    # total running time
    lookback         => 26,     # earlier snapshots searched for a file that comes back
    date_margin      => 14,     # days; only used when the previous snapshot has no manifest
    retries          => 5,      # per file
    retry_delay      => 10,     # seconds, doubled after each retry
    timeout          => 600,
);

# ################## end of defaults ####################

my @KINDS = (
    { kind => 'pdb',     dir => 'hash',         set => 'current',  key => 'pdb',
      pre => 'pdb', ext => '.ent.gz', rdir => 'data/structures/divided/pdb' },
    { kind => 'pdb-obs', dir => 'obs-hash',     set => 'obsolete', key => 'pdb',
      pre => 'pdb', ext => '.ent.gz', rdir => 'data/structures/obsolete/pdb' },
    { kind => 'cif',     dir => 'cif-hash',     set => 'current',  key => 'mmcif',
      pre => '',    ext => '.cif.gz', rdir => 'data/structures/divided/mmCIF' },
    { kind => 'cif-obs', dir => 'cif-obs-hash', set => 'obsolete', key => 'mmcif',
      pre => '',    ext => '.cif.gz', rdir => 'data/structures/obsolete/mmCIF' },
    { kind => 'xml',     dir => 'xml-hash',     set => 'current',  key => 'pdbml',
      pre => '',    ext => '.xml.gz', rdir => 'data/structures/divided/XML' },
    { kind => 'xml-obs', dir => 'xml-obs-hash', set => 'obsolete', key => 'pdbml',
      pre => '',    ext => '.xml.gz', rdir => 'data/structures/obsolete/XML' },
    { kind => 'bundle',  dir => 'bundle-hash',  set => 'current',  key => 'pdb_bundle',
      pre => '',    ext => '-pdb-bundle.tar.gz', rdir => 'compatible/pdb_bundle', perentry => 1 },
);

my @HOLDINGS = qw(current_file_holdings all_removed_entries
                  released_structures_last_modified_dates
                  obsolete_structures_last_modified_dates);

my @COLS = qw(path entry kind storage size md5 sha256 last_modified etag
              listed_date status source);
# files/REMOVED.tsv: the first ten manifest columns and the reason
my @RCOLS = (@COLS[0 .. 9], 'reason');

my $ID_RE = qr/[0-9][a-z0-9]{3}/;
my $STORAGE_RE = qr{\d{6}/files/[^/\t]+/[^/\t]+/[^/\t]+};

my ($data, $ua, $date, $stg, $deadline);
my %TMP;                # temporary files to remove on exit
my $verbose = 0;
my %c;                  # counters: $c{kind}{what}
my %n = (retries => 0, resumed => 0, repaired => 0, notkept => 0, verified => 0,
         unverified => 0, mp_unchecked => 0);
my @failed;             # [what, message]
my @notes;              # lines for the summary
my $interrupted;
my $abort;              # reason to stop the run early
my (%headCache, $headFh);
my %fetched;            # {path} = { len, etag, md5 } of files fetched into the staging directory
my %refetch;
my %syncDirs;           # directories of the staging directory this run wrote into

sub usage {
    print <<"EOT";
usage: snapshot.pl [options]
       snapshot.pl --manifest-only SNAPSHOT [--manifest-out FILE] [--reference MANIFEST ...]
       snapshot.pl --manifest-only SNAPSHOT --verify [--manifest-out FILE]
       snapshot.pl --adopt SNAPSHOT

  --root DIR             archive root (PAST_ROOT; $opt{root})
  --source URL           wwPDB base URL (PAST_SOURCE; $opt{source})
  --tmpdir DIR           scratch space (PAST_TMPDIR; $opt{tmpdir})
  --onhold-url URL       on-hold sequences (PAST_ONHOLD_URL)
  --seqres-path PATH     sequence file below the base URL ($opt{seqres_path})
  --index-path PATH      index directory below the base URL ($opt{index_path})
  --index-required LIST  comma-separated index files that must be present
  --date YYMMDD          snapshot name (PAST_DATE; today)
  --init                 first snapshot in an archive without snapshots (never
                         prompts); an interrupted first run is resumed
  --dry-run              compare holdings with the previous snapshot and stop
  --list FILE            write the names of new, candidate, removed and
                         reappeared files to FILE ("-": standard output)
  --check-all            HEAD every file, not only those whose listed date changed
  --hash-reused          compute missing checksums of reused files in this run
  --verify               read every reused file and compare it with the checksums
                         in the manifest; with --manifest-only, every file of
                         that snapshot
  --fresh                discard files left by an interrupted run
  --refetch PATH         download this file and compare it with the archived file
                         as it is on disk; a damaged archived file is replaced
                         (path as in the manifest; may be repeated)
  --etag-not-md5         do not treat any ETag as the MD5 of the file
  --adopt SNAPSHOT       make 'latest' point to a complete snapshot newer than it
                         (or to the newest one, when there is no 'latest')
 guards, checked before anything is created (per kind of file):
  --max-count-change F   change in the number of files ($opt{max_count_change})
  --max-deletions N      files that disappear ($opt{max_deletions})
  --max-candidates F     files to check with HEAD ($opt{max_candidates})
  --min-pdb-share F      current entries with a PDB-format file ($opt{min_pdb_share})
  --max-still-listed N   removed files whose entry is still listed ($opt{max_still_listed})
  --max-vanished N       removed files whose entry is listed nowhere ($opt{max_vanished})
  --max-reappeared N     new files that an earlier snapshot removed ($opt{max_reappeared})
  --allow-new-kind       accept files of a kind the previous snapshot had none of
 limits during the run:
  --max-changed F        files stored with changed content, per kind ($opt{max_changed})
  --max-repaired N       missing archived files fetched again ($opt{max_repaired})
  --max-failures N       failed files ($opt{max_failures})
  --max-compared N       files downloaded to compare and not kept, because they
                         are the archived copy ($opt{max_compared})
  --max-hours H          running time; 0 for no limit ($opt{max_hours})
 other:
  --lookback N           earlier snapshots searched for reappearing files ($opt{lookback})
  --date-margin DAYS     see README ($opt{date_margin})
  --retries N            per file ($opt{retries})
  --retry-delay SEC      first delay, doubled each retry ($opt{retry_delay})
  --timeout SEC          HTTP timeout ($opt{timeout})
  --gzip PROGRAM         ($opt{gzip})
  --verbose              one line per file
EOT
    exit($_[0] // 0);
}

sub info { print @_, "\n"; }
sub logv { print @_, "\n" if $verbose; }

END {
    my $status = $?;
    unlink(grep { -e $_ } keys %TMP) if %TMP;
    $? = $status;
}

MAIN : {
    $| = 1;
    my (@reference, @refetch);
    GetOptions(
        'root=s'             => \$opt{root},
        'source=s'           => \$opt{source},
        'tmpdir=s'           => \$opt{tmpdir},
        'onhold-url=s'       => \$opt{onhold_url},
        'seqres-path=s'      => \$opt{seqres_path},
        'index-path=s'       => \$opt{index_path},
        'index-required=s'   => \$opt{index_required},
        'date=s'             => \$opt{date},
        'init'               => \$opt{init},
        'dry-run'            => \$opt{dry_run},
        'list=s'             => \$opt{list},
        'check-all'          => \$opt{check_all},
        'hash-reused'        => \$opt{hash_reused},
        'verify'             => \$opt{verify},
        'fresh'              => \$opt{fresh},
        'refetch=s'          => \@refetch,
        'etag-not-md5'       => \$opt{etag_not_md5},
        'adopt=s'            => \$opt{adopt},
        'allow-new-kind'     => \$opt{allow_new_kind},
        'max-count-change=f' => \$opt{max_count_change},
        'max-deletions=i'    => \$opt{max_deletions},
        'max-candidates=f'   => \$opt{max_candidates},
        'min-pdb-share=f'    => \$opt{min_pdb_share},
        'max-still-listed=i' => \$opt{max_still_listed},
        'max-vanished=i'     => \$opt{max_vanished},
        'max-reappeared=i'   => \$opt{max_reappeared},
        'max-changed=f'      => \$opt{max_changed},
        'max-repaired=i'     => \$opt{max_repaired},
        'max-failures=i'     => \$opt{max_failures},
        'max-compared=i'     => \$opt{max_compared},
        'max-hours=f'        => \$opt{max_hours},
        'lookback=i'         => \$opt{lookback},
        'date-margin=i'      => \$opt{date_margin},
        'retries=i'          => \$opt{retries},
        'retry-delay=f'      => \$opt{retry_delay},
        'timeout=i'          => \$opt{timeout},
        'gzip=s'             => \$opt{gzip},
        'manifest-only=s'    => \$opt{manifest_only},
        'manifest-out=s'     => \$opt{manifest_out},
        'reference=s'        => \@reference,
        'verbose'            => \$verbose,
        'help'               => sub { usage(0) },
    ) or usage(1);
    usage(1) if @ARGV;

    $opt{source} =~ s{/+$}{};
    $data = "$opt{root}/data";
    $refetch{$_} = 1 foreach (@refetch);
    $deadline = $opt{max_hours} > 0 ? $^T + 3600 * $opt{max_hours} : undef;
    $SIG{INT} = $SIG{TERM} = $SIG{HUP} = sub {
        $interrupted = $_[0];
        die("Interrupted by signal $_[0]\n");
    };
    # output to a closed pipe must not end the run
    $SIG{PIPE} = 'IGNORE';

    if (defined $opt{manifest_only}) {
        exit(manifest_only($opt{manifest_only}, \@reference));
    }

    die("tmpdir $opt{tmpdir} is not a writeable directory\n")
        if (! -d $opt{tmpdir}) || (! -w $opt{tmpdir});

    # The lock comes first, so that the state read below is the state
    # this run works on.  A dry run writes nothing and takes no lock.  In
    # a new archive there is no data directory to hold the lock yet: it
    # is made, locked and read again once the guards have passed.
    my $locked = 0;
    if (! $opt{dry_run} && -d $data) {
        lock_archive();
        $locked = 1;
    }

    # what is in the archive
    my ($snaps, $partials) = archive_state();
    my @snaps = @$snaps;
    my @partials = @$partials;

    # previous snapshot
    my $firstRun = 0;
    my $prev;
    if (! -l "$data/latest") {
        die("$data/latest must be a symlink, not a directory.\n") if (-e "$data/latest");
        if (defined $opt{adopt}) {
            # a first run that published its snapshot but could not make 'latest'
            die("--adopt cannot be used with --init\n") if $opt{init};
            adopt($opt{adopt}, undef, \@snaps);
            exit(0);
        }
        die("$data/latest does not exist.  If this is the first snapshot in a new\n".
            "archive, run again with --init (also after a first run that failed).\n".
            "If a first run published its snapshot but could not make 'latest', use\n".
            "--adopt.  Otherwise the archive is not mounted or is damaged; nothing\n".
            "was done.\n") unless $opt{init};
        # a staging directory left by a failed first run is resumed (or
        # discarded with --fresh); a dated directory is a snapshot
        die("--init given, but $data already holds ".scalar(@snaps)." snapshot(s) ($snaps[0]".
            (@snaps > 1 ? ", ..." : "")."); nothing was done.\n".
            "If the first run published $snaps[-1] but could not make 'latest', run\n".
            "snapshot.pl --adopt $snaps[-1]\n") if @snaps;
        $firstRun = 1;
    }
    else {
        die("--init given, but $data/latest exists\n") if $opt{init};
        $prev = readlink("$data/latest");
        die("Cannot read link $data/latest: $!\n") unless defined $prev;
        $prev =~ s{/+$}{};
        die("$data/latest points to '$prev', which is not a snapshot name\n")
            unless defined snapshot_day($prev);
        die("Previous snapshot $data/$prev has no link tree\n")
            unless (-d "$data/$prev/snapshot") && (-d "$data/$prev/files");
    }

    if (defined $opt{adopt}) {
        adopt($opt{adopt}, $prev, \@snaps);
        exit(0);
    }
    die("--verify needs the manifest checksums of a previous snapshot; there is none ".
        "on a first run\n") if ($opt{verify} && $firstRun);

    $date = defined $opt{date} ? $opt{date} : strftime("%y%m%d", localtime(time()));
    die("Bad snapshot date '$date' (need YYMMDD)\n") unless defined snapshot_day($date);

    if (! $firstRun) {
        die("Snapshot $date already exists and is the latest; nothing to do.\n") if ($prev eq $date);
        die("Previous snapshot $prev is not older than $date\n") if ($prev gt $date);
        # a snapshot newer than 'latest' was finished but not switched to,
        # or 'latest' was deliberately moved back
        my @newer = grep { $_ gt $prev } @snaps;
        die("Snapshot ".join(", ", @newer)." in $data is newer than 'latest' ($prev).\n".
            "Nothing was done.  Either make it the latest (snapshot.pl --adopt $newer[-1];\n".
            "it must be complete) or move it out of $data, and run again.\n") if @newer;
    }

    clean_tmpdir() if $locked;

    $ua = LWP::UserAgent->new(timeout => $opt{timeout}, keep_alive => 4,
                              agent => "PAST-snapshot/2.0");

    # holdings
    my $tmp = tempdir("past.XXXXXXXX", DIR => $opt{tmpdir}, CLEANUP => 1);
    foreach my $h (@HOLDINGS) {
        my ($r, $msg) = with_retries("holdings $h", sub {
            fetch_once("$opt{source}/holdings/$h.json.gz", "$tmp/$h.json.gz"); });
        die("Failed to download holdings file $h: $msg\n") unless $r;
    }
    my ($new, $hstat) = parse_holdings($tmp);

    # previous membership
    my ($prevMem, $prow, $prevSource) = ({}, {}, "none (first run)");
    ($prevMem, $prow, $prevSource) = read_previous($prev) unless $firstRun;
    my $cutoff = $firstRun ? undef : day_minus($prev, $opt{date_margin});

    # plan
    my (%plan, %newPath);
    foreach my $k (@KINDS) {
        my $kind = $k->{kind};
        my $p = $plan{$kind} = { new => [], cand => [], keep => [], removed => [], why => {} };
        my $nw = $new->{$kind} || {};
        my $old = $prevMem->{$kind} || {};
        foreach my $e (sort keys %$nw) {
            if (! exists $old->{$e}) {
                push(@{$p->{new}}, $e);
                $newPath{member_path($k, $e)} = 1;
            }
            elsif (is_candidate($nw->{$e}, $prow->{member_path($k, $e)}, $cutoff)) {
                push(@{$p->{cand}}, $e);
            }
            else {
                push(@{$p->{keep}}, $e);
            }
        }
        # removed files, and why: the entry moved to the other holdings
        # file (obsoleted, usually), is still in the same one (it lost
        # this format, or the file moved), or is in neither
        my $other = ($k->{set} eq 'current') ? 'obsolete' : 'current';
        foreach my $e (sort grep { ! exists $nw->{$_} } keys %$old) {
            push(@{$p->{removed}}, $e);
            $p->{why}{$e} = $hstat->{ids}{$k->{set}}{$e} ? 'still-listed'
                          : $hstat->{ids}{$other}{$e}    ? 'moved'
                          :                                'vanished';
        }
    }
    # new files that an earlier snapshot listed as removed
    my $reap = $firstRun ? {} : find_removed_earlier(\@snaps, $prev, \%newPath);
    %newPath = ();

    info("PAST snapshot $date", ($opt{dry_run} ? " (dry run)" : ""));
    info("previous: ", ($firstRun ? "none (first run)" : "$prev, membership from $prevSource"));
    info(sprintf("%-8s %9s %9s %8s %10s %8s %8s", "kind", "previous", "listed",
                 "new", "candidate", "removed", "skipped"));
    foreach my $k (@KINDS) {
        my $kind = $k->{kind};
        my $p = $plan{$kind};
        info(sprintf("%-8s %9d %9d %8d %10d %8d %8d", $kind,
                     scalar(keys %{$prevMem->{$kind} || {}}),
                     scalar(keys %{$new->{$kind} || {}}),
                     scalar(@{$p->{new}}), scalar(@{$p->{cand}}),
                     scalar(@{$p->{removed}}), $hstat->{skipped}{$kind} || 0));
    }
    info("('candidate' files get a HEAD request; 'skipped' are listed at an unexpected path)");
    info("holdings: $hstat->{entries}{current} current entries, ".
         "$hstat->{entries}{obsolete} removed entries");
    info("note: $_") foreach (@{$hstat->{notes}});
    # the removed files that need a decision, by name
    foreach my $k (@KINDS) {
        my $p = $plan{$k->{kind}};
        foreach my $e (@{$p->{removed}}) {
            next if ($p->{why}{$e} eq 'moved');
            info("removed, entry ", ($p->{why}{$e} eq 'still-listed'
                                     ? "still listed (not written to 'removed')"
                                     : "listed nowhere"), ": ", member_path($k, $e));
        }
    }
    info("reappeared (removed in $reap->{$_}{snap}): $_") foreach (sort keys %$reap);
    # what the expensive options will read or ask for
    if (($opt{check_all} || $opt{verify}) && ! $firstRun) {
        my ($nOld, $nRec, $bytes) = (0, 0, 0);
        foreach my $k (@KINDS) {
            my $p = $plan{$k->{kind}};
            foreach my $e (@{$p->{cand}}, @{$p->{keep}}) {
                $nOld++;
                my $line = $prow->{member_path($k, $e)} or next;
                my @f = split(/\t/, $line, -1);
                next if ($f[5] eq '-' && $f[6] eq '-');
                $nRec++;
                $bytes += $f[4] if ($f[4] ne '-');
            }
        }
        info("--check-all: a HEAD request for each of the $nOld files already archived.  A file ".
             "that the answer does not identify is downloaded and compared; the run stops when ".
             "more than $opt{max_compared} such downloads were not needed (--max-compared) or ".
             "more files changed than --max-changed allows.") if $opt{check_all};
        info("--verify: $nRec archived files ($bytes bytes) will be read and compared with the ".
             "checksums in the previous manifest; ".($nOld - $nRec)." have no recorded checksum ".
             "and cannot be verified.") if $opt{verify};
    }

    # guards
    my @guards = check_guards($new, $prevMem, \%plan, $hstat, $firstRun, $reap);
    my $guardsOk = 1;
    foreach my $g (@guards) {
        info(sprintf("guard %-4s %s", ($g->[0] ? "ok" : "FAIL"), $g->[1]));
        $guardsOk = 0 unless $g->[0];
    }
    write_list(\%plan, $reap) if defined $opt{list};
    if ($opt{dry_run}) {
        info("Dry run: nothing fetched, nothing created.");
        exit($guardsOk ? 0 : 2);
    }
    die("Guard failed; nothing was created and 'latest' is unchanged.\n") unless $guardsOk;

    if (! $locked) {
        # first run in a new archive: there was no data directory to lock
        make_dir($data);
        lock_archive();
        die("$data/latest appeared meanwhile; nothing was done.\n")
            if (-e "$data/latest" || -l "$data/latest");
        (my $now, $partials) = archive_state();
        die("$data holds a snapshot ($now->[0]) that appeared meanwhile; nothing was done.\n")
            if @$now;
        @partials = @$partials;
        clean_tmpdir();
    }

    # staging directory
    $stg = "$data/.$date.partial";
    prepare_staging(\@partials);
    make_dir("$stg/files/holdings");
    foreach my $h (@HOLDINGS) {
        copy("$tmp/$h.json.gz", "$stg/files/holdings/$h.json.gz")
            or die("Cannot copy $h.json.gz into the snapshot: $!\n");
        fsync_path("$stg/files/holdings/$h.json.gz");
    }
    load_head_cache();

    # sequences and index files first: they are small, and a failure
    # there should not come after all the other work
    my @rows;       # manifest lines
    my @aux;
    push(@aux, get_seqres());
    push(@aux, get_onhold()) unless @failed;
    push(@aux, make_pdbprot()) unless @failed;
    push(@aux, get_index()) unless @failed;

    # structure files
    my %staged;      # "dir/hh/name" of files that belong to this run
    my @reappeared;  # [path, removed in, outcome]
    if (! @failed) {
        KIND: foreach my $k (@KINDS) {
            my $kind = $k->{kind};
            my $p = $plan{$kind};
            my %state;
            $state{$_} = 'new'  foreach (@{$p->{new}});
            $state{$_} = 'cand' foreach (@{$p->{cand}});
            $state{$_} = 'keep' foreach (@{$p->{keep}});
            my $maxChanged = $opt{max_changed} * scalar(keys %{$prevMem->{$kind} || {}});
            # in the order of the manifest
            foreach my $e (sort { member_path($k, $a) cmp member_path($k, $b) } keys %state) {
                my $path = member_path($k, $e);
                my $row = process_entry($k, $e, $state{$e}, $new->{$kind}{$e}, $prev,
                                        prev_row($prow, $path), $reap->{$path}, \@reappeared);
                if ($row) {
                    push(@rows, join("\t", map { defined $row->{$_} ? $row->{$_} : '-' } @COLS));
                    $staged{$path} = 1 if ($row->{status} eq 'fetched');
                }
                if (! $firstRun && ($c{$kind}{changed} || 0) > $maxChanged) {
                    $abort = sprintf("more than %d %s files have changed content (--max-changed %s)",
                                     $maxChanged, $kind, $opt{max_changed});
                }
                $abort = "more than $opt{max_repaired} archived files are missing (--max-repaired)"
                    if ($n{repaired} > $opt{max_repaired});
                $abort = "more than $opt{max_compared} files were downloaded only to find that ".
                    "they are the archived copy (--max-compared).  If the server no longer ".
                    "gives usable ETags, the whole archive would be downloaded this way"
                    if ($n{notkept} > $opt{max_compared});
                $abort = "more than $opt{max_failures} files failed (--max-failures)"
                    if (@failed > $opt{max_failures});
                $abort = "running longer than $opt{max_hours} hours (--max-hours)"
                    if (defined $deadline && time() > $deadline);
                last KIND if $abort;
            }
            $c{$kind}{removed} = scalar(@{$p->{removed}});
        }
    }

    if (@failed || $abort) {
        summary();
        info("STOPPED: $abort") if $abort;
        info("FAILED: $date was not published; 'latest' is unchanged.");
        info("Files already fetched are kept in $stg and are used by the next run.");
        exit(1);
    }
    push(@rows, map { join("\t", map { defined $_ ? $_ : '-' } @$_{@COLS}) } @aux);

    # finish the staging directory
    sweep_staging(\%staged);
    write_removed(\%plan, $prev, $prow);
    if (@reappeared) {
        my $f = "$stg/files/REAPPEARED.tsv";
        open(my $fh, '>', $f) or die("Cannot write $f: $!\n");
        print $fh "#path\tremoved_in\toutcome\n";
        print $fh join("\t", @$_), "\n" foreach (@reappeared);
        close_synced($fh, $f);
    }
    foreach my $l ('seq', 'index') {
        symlink("../files/$l", "$stg/snapshot/$l")
            or die("Cannot make link $stg/snapshot/$l: $!\n");
    }

    # verify before publishing
    my $want = scalar(@aux);
    $want += scalar(keys %{$new->{$_->{kind}} || {}}) foreach (@KINDS);
    die("Internal check failed: ".scalar(@rows)." manifest rows for $want files\n")
        if (@rows != $want);
    my $pending = 0;
    foreach my $line (@rows) {
        my @f = split(/\t/, $line, -1);
        $pending++ if ($f[5] eq '-' || $f[6] eq '-');
        next unless ($f[10] eq 'fetched');
        my $size = -s "$stg/files/".substr($f[3], 13);
        die("Internal check failed: $f[3] is missing or has the wrong size\n")
            unless (defined $size && $size == $f[4]);
    }

    write_manifest("$stg/files/MANIFEST.tsv", $date, $prev, \@rows);
    close($headFh) if $headFh;
    unlink("$stg/HEADS.tsv") or $!{ENOENT} or die("Cannot remove $stg/HEADS.tsv: $!\n");

    my $marker = "$stg/COMPLETE";
    open(my $fh, '>', $marker) or die("Cannot write $marker: $!\n");
    print $fh "snapshot\t$date\n", "previous\t", ($prev // '-'), "\n",
        "finished\t", iso_time(time()), "\n", "files\t", scalar(@rows), "\n";
    close_synced($fh, $marker);
    # every directory this run wrote into, deepest first
    $syncDirs{$_} = 1 foreach ("$stg/files", "$stg/snapshot", $stg);
    fsync_dir($_) foreach (sort { length($b) <=> length($a) || $a cmp $b } keys %syncDirs);

    # publish
    rename($stg, "$data/$date") or die("Cannot rename $stg to $data/$date: $!\n");
    fsync_dir($data);
    switch_latest($date);
    root_links($firstRun);

    summary();
    info("manifest: $data/$date/files/MANIFEST.tsv, ", scalar(@rows), " rows",
         ($pending ? ", $pending without checksums (run: snapshot.pl --manifest-only $date)" : ""));
    info("OK: snapshot $date published; 'latest' -> $date");
    exit(0);
}

# ---------------------------------------------------------------------
# names and dates

sub entry_name {
    my ($k, $e) = @_;
    return $k->{pre}.$e.$k->{ext};
}

sub member_path {
    my ($k, $e) = @_;
    return "$k->{dir}/".substr($e, 1, 2)."/".entry_name($k, $e);
}

sub remote_path {
    my ($k, $e) = @_;
    my $hh = substr($e, 1, 2);
    return $k->{perentry} ? "$k->{rdir}/$hh/$e/".entry_name($k, $e)
                          : "$k->{rdir}/$hh/".entry_name($k, $e);
}

sub all_dir {
    my ($k) = @_;
    (my $a = $k->{dir}) =~ s/hash$/all/;
    return $a;
}

sub name_entry {
    my ($k, $name) = @_;
    return ($name =~ /^\Q$k->{pre}\E($ID_RE)\Q$k->{ext}\E$/) ? $1 : undef;
}

# YYMMDD -> epoch at 00:00 UTC, or undef
sub snapshot_day {
    my ($s) = @_;
    return undef unless (defined $s && $s =~ /^(\d\d)(\d\d)(\d\d)$/);
    my $t = eval { timegm(0, 0, 0, $3, $2 - 1, 2000 + $1) };
    return $t;
}

sub day_minus {
    my ($snap, $days) = @_;
    return strftime("%Y%m%d", gmtime(snapshot_day($snap) - 86400 * $days));
}

# listed date -> YYYYMMDD; dies when it cannot be parsed
sub day_key {
    my ($d, $where) = @_;
    if (defined $d && !ref($d) &&
        $d =~ /^(\d{4})-(\d\d)-(\d\d)(?:T\d\d:\d\d:\d\d(?:\.\d+)?(?:Z|[+-]\d\d:?\d\d)?)?$/ &&
        $2 >= 1 && $2 <= 12 && $3 >= 1 && $3 <= 31) {
        return "$1$2$3";
    }
    die("Unparseable date '".(defined $d ? $d : 'undef')."' $where; nothing was done.\n");
}

sub iso_time {
    return strftime("%Y-%m-%dT%H:%M:%SZ", gmtime($_[0]));
}

sub iso_epoch {
    my ($s) = @_;
    return undef unless (defined $s && $s =~ /^(\d{4})-(\d\d)-(\d\d)T(\d\d):(\d\d):(\d\d)Z$/);
    return eval { timegm($6, $5, $4, $3, $2 - 1, $1) };
}

# ---------------------------------------------------------------------
# archive state

sub make_dir {
    my ($d) = @_;
    return if (-d $d);
    make_path($d, { mode => 0755, error => \my $err });
    die("Cannot create directory $d\n") if (@$err || ! -d $d);
}

sub fsync_path {
    my ($f) = @_;
    open(my $fh, '<', $f) or die("Cannot open $f: $!\n");
    $fh->sync or die("Cannot sync $f: $!\n");
    close($fh);
}

# best effort: not every file system syncs a directory
sub fsync_dir {
    my ($d) = @_;
    if (open(my $fh, '<', $d)) {
        $fh->sync;
        close($fh);
    }
}

sub close_synced {
    my ($fh, $f) = @_;
    ($fh->flush && $fh->sync && close($fh)) or die("Cannot write $f: $!\n");
}

# Dated directories (snapshots) and staging directories in the data
# directory, by name.
sub archive_state {
    my (@snaps, @partials);
    if (-d $data) {
        opendir(my $dh, $data) or die("Cannot read $data: $!\n");
        foreach my $name (readdir($dh)) {
            push(@snaps, $name) if ($name =~ /^\d{6}$/ && -d "$data/$name" && ! -l "$data/$name");
            push(@partials, $name) if ($name =~ /^\.\d{6}\.partial$/);
        }
        closedir($dh);
    }
    return ([ sort @snaps ], [ sort @partials ]);
}

my $lockFh;
sub lock_archive {
    return if $lockFh;
    my $f = "$data/.snapshot.lock";
    open($lockFh, '>>', $f) or die("Cannot open lock file $f: $!\n");
    flock($lockFh, LOCK_EX | LOCK_NB)
        or die("Another snapshot.pl holds the lock $f; nothing was done.\n");
}

# Temporary directories of runs that were killed.  Only called while
# holding the lock, so no run on this archive is using one.
sub clean_tmpdir {
    opendir(my $dh, $opt{tmpdir}) or return;
    my @old = grep { /^past\.[A-Za-z0-9_]{8}$/ } readdir($dh);
    closedir($dh);
    foreach my $name (@old) {
        my $d = "$opt{tmpdir}/$name";
        next unless (-d $d && ! -l $d && -o $d && -M $d > 1);
        remove_tree($d, { safe => 1 });
    }
}

# only called while holding the lock
sub switch_latest {
    my ($to) = @_;
    my $tmp = "$data/latest.tmp";
    unlink($tmp);
    symlink($to, $tmp) or die("Cannot make link $tmp: $!\n");
    if (! rename($tmp, "$data/latest")) {
        my $err = $!;
        unlink($tmp);
        die("Cannot switch $data/latest to $to: $err\n");
    }
    fsync_dir($data);
}

# --adopt: make a finished snapshot the latest.  $prev: what 'latest'
# points to; undef when there is no 'latest' (a first run published its
# snapshot and could not make the link), and then only the newest
# snapshot can be adopted.
sub adopt {
    my ($snap, $prev, $snaps) = @_;
    die("'$snap' is not a snapshot name (YYMMDD)\n") unless defined snapshot_day($snap);
    if (defined $prev) {
        die("$snap is not newer than the latest snapshot ($prev)\n") unless ($snap gt $prev);
    }
    else {
        my @newer = grep { $_ gt $snap } @$snaps;
        die("There is no 'latest', and $snap is not the newest snapshot in $data (".
            join(", ", @newer)." is newer); nothing was done.\n") if @newer;
    }
    my $d = "$data/$snap";
    die("$d is not a directory\n") if (-l $d || ! -d $d);
    die("$d has no COMPLETE marker: it was not made by a finished run\n") unless (-f "$d/COMPLETE");
    foreach my $k (@KINDS) {
        foreach my $sub ($k->{dir}, all_dir($k)) {
            die("$d has no link tree ($sub)\n") unless (-d "$d/snapshot/$sub");
        }
    }
    my $rows = read_manifest("$d/files/MANIFEST.tsv");
    lock_archive();
    switch_latest($snap);
    info("'latest' -> $snap (", scalar(keys %$rows), " files in its manifest)");
    # the run that published it did not come to this
    root_links(0);
    info("note: $_") foreach (@notes);
}

sub root_links {
    my ($firstRun) = @_;
    foreach my $l ('seq', 'index', map { ($_->{dir}, all_dir($_)) } @KINDS) {
        next if (-l "$opt{root}/$l");
        if (! symlink("data/latest/snapshot/$l", "$opt{root}/$l")) {
            my $msg = "Cannot make link $opt{root}/$l: $!";
            die("$msg\n") if $firstRun;
            push(@notes, $msg);
        }
    }
}

sub prepare_staging {
    my ($partials) = @_;
    my @stale = map { "$data/$_" } sort @$partials;
    if ($opt{fresh}) {
        foreach my $d (@stale) {
            remove_tree($d, { safe => 1 });
            die("Cannot remove $d\n") if (-e $d);
            info("removed $d (--fresh)");
        }
        @stale = ();
    }
    die("More than one interrupted run found (".join(", ", @stale).
        "); remove all but one, or use --fresh.\n") if (@stale > 1);
    if (@stale) {
        if ($stale[0] ne $stg) {
            rename($stale[0], $stg) or die("Cannot rename $stale[0] to $stg: $!\n");
        }
        info("resuming: files already fetched into ", basename($stale[0]),
             " are checked and kept if they match");
        # everything except the fetched structure files is made again
        foreach my $d ("$stg/snapshot", "$stg/files/seq", "$stg/files/index",
                       "$stg/files/holdings") {
            next unless (-e $d);
            remove_tree($d, { safe => 1 });
            die("Cannot remove $d\n") if (-e $d);
        }
        # a list left here would be published with the new snapshot
        foreach my $f ("$stg/files/MANIFEST.tsv", "$stg/files/REMOVED.tsv",
                       "$stg/files/REAPPEARED.tsv", "$stg/COMPLETE") {
            unlink($f) or $!{ENOENT} or die("Cannot remove $f: $!\n");
        }
    }
    my @dirs = ("$stg/files/index", "$stg/files/seq", "$stg/files/holdings");
    foreach my $k (@KINDS) {
        push(@dirs, "$stg/files/$k->{dir}", "$stg/snapshot/$k->{dir}", "$stg/snapshot/".all_dir($k));
    }
    foreach my $d (@dirs) {
        make_dir($d);
        $syncDirs{$d} = 1;
    }
}

# Remove everything under the staged files/<kind> directories that this
# run did not fetch: files of an interrupted run that are no longer
# wanted, and temporary files.  ParsePDB reads all of files/hash.
sub sweep_staging {
    my ($staged) = @_;
    my $removed = 0;
    foreach my $k (@KINDS) {
        my $top = "$stg/files/$k->{dir}";
        opendir(my $dh, $top) or die("Cannot read $top: $!\n");
        my @subs = grep { $_ ne '.' && $_ ne '..' } readdir($dh);
        closedir($dh);
        foreach my $s (@subs) {
            if (-d "$top/$s" && ! -l "$top/$s") {
                opendir(my $sh, "$top/$s") or die("Cannot read $top/$s: $!\n");
                my @files = grep { $_ ne '.' && $_ ne '..' } readdir($sh);
                closedir($sh);
                my $left = 0;
                foreach my $f (@files) {
                    if ($staged->{"$k->{dir}/$s/$f"}) {
                        $left++;
                        next;
                    }
                    unlink("$top/$s/$f") or die("Cannot remove $top/$s/$f: $!\n");
                    $removed++;
                    $syncDirs{"$top/$s"} = 1;
                }
                if (! $left && rmdir("$top/$s")) {
                    delete $syncDirs{"$top/$s"};
                }
            }
            else {
                unlink("$top/$s") or die("Cannot remove $top/$s: $!\n");
                $removed++;
            }
        }
    }
    push(@notes, "$removed stale file(s) of an interrupted run removed") if $removed;
}

# files/<dir>/removed, as before: the names ParsePDB marks obsolete (for
# hash and obs-hash).  A file whose entry is still listed in the same
# holdings file is not written there.  files/REMOVED.tsv records every
# removed file with where its content is kept, so that a file that
# comes back can be recognised.
sub write_removed {
    my ($plan, $prev, $prow) = @_;
    my @lines;
    foreach my $k (@KINDS) {
        my $p = $plan->{$k->{kind}};
        my @gone = @{$p->{removed}};
        next unless @gone;
        my @listed = grep { $p->{why}{$_} ne 'still-listed' } @gone;
        if (@listed) {
            my $f = "$stg/files/$k->{dir}/removed";
            open(my $fh, '>', $f) or die("Cannot write $f: $!\n");
            print $fh entry_name($k, $_), "\n" foreach (@listed);
            close_synced($fh, $f);
        }
        foreach my $e (@gone) {
            my $path = member_path($k, $e);
            my $r = prev_row($prow, $path);
            if (! $r) {
                my $l = readlink("$data/$prev/snapshot/$path");
                $r = { path => $path, entry => $e, kind => $k->{kind},
                       storage => ((defined $l && $l =~ m{^(?:\.\./){4}($STORAGE_RE)$}) ? $1 : '-') };
            }
            push(@lines, join("\t", (map { defined $r->{$_} ? $r->{$_} : '-' } @RCOLS[0 .. 9]),
                              $p->{why}{$e}));
            push(@notes, "removed but not written to '$k->{dir}/removed' (entry still listed): ".
                 entry_name($k, $e)) if ($p->{why}{$e} eq 'still-listed');
        }
    }
    return unless @lines;
    my $f = "$stg/files/REMOVED.tsv";
    open(my $fh, '>', $f) or die("Cannot write $f: $!\n");
    print $fh "#", join("\t", @RCOLS), "\n";
    print $fh "$_\n" foreach (@lines);
    close_synced($fh, $f);
}

# New files that one of the last --lookback snapshots listed as removed.
# Returns {path}{snap, row}; row->{storage} is '-' when the kept copy
# cannot be located (snapshots made before REMOVED.tsv existed, with a
# pruned predecessor).
sub find_removed_earlier {
    my ($snaps, $prev, $wanted) = @_;
    my %found;
    return \%found unless %$wanted;
    my @all = grep { $_ le $prev } @$snaps;
    my $first = @all - $opt{lookback};
    $first = 0 if ($first < 0);
    for (my $i = $#all; $i >= $first; $i--) {
        my $s = $all[$i];
        my $f = "$data/$s/files/REMOVED.tsv";
        if (open(my $fh, '<', $f)) {
            while (my $line = <$fh>) {
                chomp($line);
                next if ($line =~ /^#/ || $line eq '');
                my @v = split(/\t/, $line, -1);
                die("$f line $.: expected ".scalar(@RCOLS)." columns\n") if (@v != @RCOLS);
                # as for a manifest row; "-" when the kept copy was not found
                die("$f line $.: unexpected storage path for $v[0]\n")
                    unless ($v[3] eq '-' || $v[3] =~ /^$STORAGE_RE$/);
                next unless ($wanted->{$v[0]} && ! $found{$v[0]});
                my %r;
                @r{@RCOLS} = @v;
                $found{$v[0]} = { snap => $s, row => \%r };
            }
            close($fh);
            next;
        }
        die("Cannot read $f: $!\n") unless $!{ENOENT};
        foreach my $k (@KINDS) {
            $f = "$data/$s/files/$k->{dir}/removed";
            if (! open(my $fh, '<', $f)) {
                die("Cannot read $f: $!\n") unless $!{ENOENT};
                next;
            }
            else {
                while (my $name = <$fh>) {
                    chomp($name);
                    my $e = name_entry($k, $name);
                    next unless defined $e;
                    my $path = member_path($k, $e);
                    next unless ($wanted->{$path} && ! $found{$path});
                    my $storage = '-';
                    if ($i > 0) {
                        my $l = readlink("$data/$all[$i - 1]/snapshot/$path");
                        $storage = $1 if (defined $l && $l =~ m{^(?:\.\./){4}($STORAGE_RE)$});
                    }
                    $found{$path} = { snap => $s, row => {
                        path => $path, entry => $e, kind => $k->{kind}, storage => $storage,
                        size => '-', md5 => '-', sha256 => '-', last_modified => '-', etag => '-' } };
                }
                close($fh);
            }
        }
    }
    return \%found;
}

# ---------------------------------------------------------------------
# holdings

sub load_json {
    my ($file) = @_;
    my $fh = open_cmd($opt{gzip}, '-dc', '--', $file);
    my $text;
    {
        local $/;
        $text = <$fh>;
    }
    close($fh) or die("Cannot decompress $file\n");
    my $decode;
    foreach my $m ('Cpanel::JSON::XS', 'JSON::XS', 'JSON::MaybeXS', 'JSON::PP') {
        (my $f = "$m.pm") =~ s{::}{/}g;
        next unless eval { require $f; 1 };
        $decode = $m->can('decode_json');
        last if $decode;
    }
    die("No JSON module found\n") unless $decode;
    my $d = eval { $decode->($text) };
    die("Cannot parse $file: $@\n") unless (ref($d) eq 'HASH');
    return $d;
}

# Returns ({kind}{entry} = listed date or ''), statistics.
sub parse_holdings {
    my ($dir) = @_;
    my (%new, %stat);
    $stat{notes} = [];
    foreach my $set ('current', 'obsolete') {
        my $hold  = load_json("$dir/".($set eq 'current' ? $HOLDINGS[0] : $HOLDINGS[1]).".json.gz");
        my $dates = load_json("$dir/".($set eq 'current' ? $HOLDINGS[2] : $HOLDINGS[3]).".json.gz");
        day_key($dates->{$_}, "for $_ in the $set dates file") foreach (keys %$dates);
        $stat{entries}{$set} = 0;
        $stat{ids}{$set} = {};
        my ($badId, $noDate, $obsBundle) = (0, 0, 0);
        foreach my $code (keys %$hold) {
            my $e = lc $code;
            if ($e !~ /^$ID_RE$/) {
                $badId++;
                next;
            }
            my $rec = $hold->{$code};
            die("Unexpected holdings record for $code\n") unless (ref($rec) eq 'HASH');
            my $ct = ($set eq 'current') ? $rec : ($rec->{content_type} || {});
            die("Unexpected holdings record for $code\n") unless (ref($ct) eq 'HASH');
            $stat{entries}{$set}++;
            $stat{ids}{$set}{$e} = 1;
            $obsBundle++ if ($set eq 'obsolete' && exists $ct->{pdb_bundle});
            my $listed = defined $dates->{$code} ? $dates->{$code} : '';
            my $any = 0;
            foreach my $k (grep { $_->{set} eq $set } @KINDS) {
                my $list = $ct->{$k->{key}};
                next unless defined $list;
                die("Unexpected holdings record for $code ($k->{key})\n")
                    unless (ref($list) eq 'ARRAY');
                my $want = "/pdb/".remote_path($k, $e);
                if (grep { $_ eq $want } @$list) {
                    $new{$k->{kind}}{$e} = $listed;
                    $any = 1;
                }
                else {
                    $stat{skipped}{$k->{kind}}++;
                }
            }
            $noDate++ if ($any && $listed eq '');
        }
        push(@{$stat{notes}}, "$badId $set entries have an id that is not four characters; not archived") if $badId;
        push(@{$stat{notes}}, "$noDate $set entries have no listed date; their files are always checked") if $noDate;
        push(@{$stat{notes}}, "$obsBundle removed entries have a bundle; obsolete bundles are not archived") if $obsBundle;
    }
    return (\%new, \%stat);
}

# ---------------------------------------------------------------------
# manifest

# Returns {path} = line.  Lines are kept as text to save memory.
sub read_manifest {
    my ($file, $partial) = @_;
    my (%rows, $trailer);
    open(my $fh, '<', $file) or die("Cannot read $file: $!\n");
    while (my $line = <$fh>) {
        my $complete = ($line =~ s/\n$//);
        if ($line =~ /^#/) {
            $trailer = $1 if ($line =~ /^#rows\t(\d+)$/);
            next;
        }
        next if ($line eq '');
        my $nf = ($line =~ tr/\t//) + 1;
        if ($nf != @COLS || ! $complete) {
            next if $partial;
            die("$file line $.: expected ".scalar(@COLS)." columns\n");
        }
        $rows{substr($line, 0, index($line, "\t"))} = $line;
    }
    close($fh);
    if (! $partial) {
        die("$file is truncated (no row count at the end)\n") unless defined $trailer;
        die("$file has ".scalar(keys %rows)." rows, but says $trailer\n")
            if ($trailer != scalar(keys %rows));
    }
    return \%rows;
}

sub row_hash {
    my ($line) = @_;
    my %r;
    @r{@COLS} = split(/\t/, $line, -1);
    return \%r;
}

sub prev_row {
    my ($prow, $path) = @_;
    my $line = $prow->{$path};
    return defined $line ? row_hash($line) : undef;
}

# $rows: lines, already in the order wanted
sub write_manifest {
    my ($file, $snap, $prev, $rows) = @_;
    my $tmp = "$file.tmp.$$";
    $TMP{$tmp} = 1;
    open(my $fh, '>', $tmp) or die("Cannot write $tmp: $!\n");
    print $fh "#PAST-manifest\t1\n", "#snapshot\t$snap\n", "#previous\t", ($prev // '-'), "\n",
        "#", join("\t", @COLS), "\n";
    foreach my $line (@$rows) {
        print $fh $line, "\n" or die("Cannot write $tmp: $!\n");
    }
    print $fh "#rows\t", scalar(@$rows), "\n";
    close_synced($fh, $tmp);
    rename($tmp, $file) or die("Cannot rename $tmp to $file: $!\n");
    delete $TMP{$tmp};
}

# Membership of the previous snapshot: from its manifest when it has one,
# otherwise from one directory read per kind (the "all" directories).
sub read_previous {
    my ($prev) = @_;
    my (%mem, $prow, $source);
    my $mf = "$data/$prev/files/MANIFEST.tsv";
    if (-f $mf) {
        $prow = read_manifest($mf);
        $source = "its manifest";
        my %byKind = map { $_->{kind} => $_ } @KINDS;
        foreach my $line (values %$prow) {
            my @f = split(/\t/, $line, -1);
            my $k = $byKind{$f[2]} or next;
            my $e = name_entry($k, basename($f[0]));
            die("$mf: unexpected path $f[0]\n")
                unless (defined $e && $e eq $f[1] && $f[0] eq member_path($k, $e));
            die("$mf: unexpected storage path for $f[0]\n") unless ($f[3] =~ /^$STORAGE_RE$/);
            day_key($f[9], "for $f[0] in $mf") if ($f[9] ne '-');
            $mem{$k->{kind}}{$e} = 1;
        }
    }
    else {
        $prow = {};
        $source = "its link tree";
        my $odd = 0;
        foreach my $k (@KINDS) {
            my $d = "$data/$prev/snapshot/".all_dir($k);
            opendir(my $dh, $d) or die("Cannot read previous membership $d: $!\n");
            while (defined(my $name = readdir($dh))) {
                next if ($name eq '.' || $name eq '..');
                my $e = name_entry($k, $name);
                if (defined $e) {
                    $mem{$k->{kind}}{$e} = 1;
                }
                else {
                    $odd++;
                }
            }
            closedir($dh);
        }
        push(@notes, "$odd name(s) in the previous link tree were not recognised and are dropped") if $odd;
    }
    my $total = 0;
    $total += scalar(keys %{$mem{$_}}) foreach (keys %mem);
    die("Previous snapshot $prev has no files in its membership; nothing was done.\n") unless $total;
    return (\%mem, $prow, $source);
}

# ---------------------------------------------------------------------
# plan and guards

# A candidate gets a HEAD request; the content identity decides.
# $pline: line of the previous manifest, if any.
sub is_candidate {
    my ($listed, $pline, $cutoff) = @_;
    return 1 if $opt{check_all};
    return 1 if ($listed eq '');
    if (defined $pline) {
        my $was = (split(/\t/, $pline, -1))[9];
        return ($was ne $listed) ? 1 : 0 if ($was ne '-');
    }
    # no record of what was listed last time: check everything listed as
    # modified since shortly before the previous snapshot
    return (day_key($listed, "in the holdings") ge $cutoff) ? 1 : 0;
}

sub check_guards {
    my ($new, $prevMem, $plan, $hstat, $firstRun, $reap) = @_;
    my @g;
    my $nCur = $hstat->{entries}{current};
    push(@g, [ $nCur > 0, "current holdings list $nCur entries" ]);
    my $nPdb = scalar(keys %{$new->{pdb} || {}});
    my $share = $nCur ? $nPdb / $nCur : 0;
    push(@g, [ $share >= $opt{min_pdb_share},
               sprintf("PDB-format share %.4f (%d of %d), floor %s",
                       $share, $nPdb, $nCur, $opt{min_pdb_share}) ]);
    return @g if $firstRun;
    my ($still, $vanished) = (0, 0);
    foreach my $k (@KINDS) {
        my $kind = $k->{kind};
        my $p = $plan->{$kind};
        my $old = scalar(keys %{$prevMem->{$kind} || {}});
        my $now = scalar(keys %{$new->{$kind} || {}});
        if (! $old) {
            push(@g, [ ($now == 0 || $opt{allow_new_kind}),
                       "$kind count 0 -> $now (the previous snapshot had none; --allow-new-kind accepts)" ]);
        }
        else {
            my $change = abs($now - $old) / $old;
            push(@g, [ $change <= $opt{max_count_change},
                       sprintf("%s count %d -> %d (change %.4f, limit %s)",
                               $kind, $old, $now, $change, $opt{max_count_change}) ]);
        }
        my $gone = scalar(@{$p->{removed}});
        push(@g, [ $gone <= $opt{max_deletions},
                   sprintf("%s removed %d (limit %d)", $kind, $gone, $opt{max_deletions}) ]);
        if (! $opt{check_all}) {
            my $cand = scalar(@{$p->{cand}});
            my $frac = $now ? $cand / $now : 0;
            push(@g, [ $frac <= $opt{max_candidates},
                       sprintf("%s candidates %d (%.4f of listed, limit %s)",
                               $kind, $cand, $frac, $opt{max_candidates}) ]);
        }
        $still    += scalar(grep { $_ eq 'still-listed' } values %{$p->{why}});
        $vanished += scalar(grep { $_ eq 'vanished' } values %{$p->{why}});
    }
    push(@g, [ $still <= $opt{max_still_listed},
               "removed files whose entry is still listed: $still (limit $opt{max_still_listed})" ]);
    push(@g, [ $vanished <= $opt{max_vanished},
               "removed files whose entry is listed nowhere: $vanished (limit $opt{max_vanished})" ]);
    my $nReap = scalar(keys %$reap);
    push(@g, [ $nReap <= $opt{max_reappeared},
               "new files that an earlier snapshot removed: $nReap (limit $opt{max_reappeared})" ]);
    return @g;
}

# --list: every name, one per line
sub write_list {
    my ($plan, $reap) = @_;
    my $fh;
    if ($opt{list} eq '-') {
        $fh = \*STDOUT;
    }
    else {
        open($fh, '>', $opt{list}) or die("Cannot write $opt{list}: $!\n");
    }
    print $fh "#action\tpath\tdetail\n";
    foreach my $k (@KINDS) {
        my $p = $plan->{$k->{kind}};
        foreach my $e (@{$p->{new}}) {
            my $path = member_path($k, $e);
            print $fh ($reap->{$path} ? "reappeared\t$path\tremoved in $reap->{$path}{snap}\n"
                                      : "new\t$path\t\n");
        }
        print $fh "candidate\t", member_path($k, $_), "\t\n" foreach (@{$p->{cand}});
        foreach my $e (@{$p->{removed}}) {
            my $why = $p->{why}{$e};
            print $fh "removed\t", member_path($k, $e), "\t",
                ($why eq 'moved' ? "entry is now in the other holdings file; written to $k->{dir}/removed"
                 : $why eq 'vanished' ? "entry is listed nowhere; written to $k->{dir}/removed"
                 : "entry is still listed; not written to $k->{dir}/removed"), "\n";
        }
    }
    if ($opt{list} ne '-') {
        close($fh) or die("Cannot write $opt{list}: $!\n");
    }
}

# ---------------------------------------------------------------------
# checksums and child processes

# Read the output of a command (list form, no shell, stderr discarded).
sub open_cmd {
    my (@cmd) = @_;
    my $pid = open(my $fh, '-|');
    die("Cannot fork: $!\n") unless defined $pid;
    if (! $pid) {
        open(STDERR, '>', '/dev/null');
        no warnings 'exec';
        exec { $cmd[0] } @cmd;
        POSIX::_exit(127);
    }
    binmode($fh);
    return $fh;
}

sub md5_file {
    my ($file) = @_;
    open(my $fh, '<', $file) or die("Cannot read $file: $!\n");
    binmode($fh);
    my $d = Digest::MD5->new->addfile($fh)->hexdigest;
    close($fh);
    return $d;
}

sub sha256_file {
    my ($file) = @_;
    open(my $fh, '<', $file) or die("Cannot read $file: $!\n");
    binmode($fh);
    my $d = Digest::SHA->new(256)->addfile($fh)->hexdigest;
    close($fh);
    return $d;
}

# SHA-256 of the decompressed content; also the gzip integrity test.
# Returns undef when gzip reports any problem.
sub gz_sha256 {
    my ($file) = @_;
    my $fh = open_cmd($opt{gzip}, '-dc', '--', $file);
    my $d = Digest::SHA->new(256)->addfile($fh)->hexdigest;
    return undef unless close($fh);
    return $d;
}

# $name decides whether the file is gzip (a temporary name does not say)
sub content_sha256 {
    my ($file, $name) = @_;
    $name = $file unless defined $name;
    return ($name =~ /\.gz$/) ? gz_sha256($file) : sha256_file($file);
}

# ---------------------------------------------------------------------
# HTTP

# Strong ETag without its quotes; undef for none, empty or weak (W/"...")
sub clean_etag {
    my ($e) = @_;
    return undef unless defined $e;
    $e =~ s/^\s+|\s+$//g;
    return undef if ($e =~ m{^W/});
    $e = $1 if ($e =~ /^"(.*)"$/);
    return undef if ($e eq '' || $e =~ /[\t\n"]/);
    return $e;
}

# the MD5 when the ETag is one (single-part object), else undef
sub etag_md5 {
    my ($e) = @_;
    return undef if $opt{etag_not_md5};
    return (defined $e && $e =~ /^[0-9a-fA-F]{32}$/) ? lc($e) : undef;
}

sub etag_multipart {
    my ($e) = @_;
    return (defined $e && $e =~ /^[0-9a-f]{32}-\d+$/) ? 1 : 0;
}

# The ETag of a multipart upload is the MD5 of the MD5s of its parts,
# "-", and the number of parts.  True when the file gives this ETag for
# one of the usual part sizes (files.wwpdb.org uses 8 MiB); then the
# content is the same.  False proves nothing.
sub multipart_matches {
    my ($file, $etag) = @_;
    return 0 if $opt{etag_not_md5};
    my ($n) = ($etag =~ /-(\d+)$/);
    my $size = -s $file;
    return 0 unless ($size && $n);
    foreach my $mb (8, 16, 5, 32, 64) {
        my $ps = $mb << 20;
        next unless (int(($size + $ps - 1) / $ps) == $n);
        open(my $fh, '<', $file) or die("Cannot read $file: $!\n");
        binmode($fh);
        my ($all, $buf) = ('');
        while (read($fh, $buf, $ps)) {
            $all .= Digest::MD5::md5($buf);
        }
        close($fh);
        return 1 if (Digest::MD5::md5_hex($all)."-$n" eq $etag);
    }
    return 0;
}

sub clean_length {
    my ($l) = @_;
    return (defined $l && $l =~ /^\d+$/) ? $l : undef;
}

# Runs $code until it returns a defined value; bounded, delay starts
# afresh for every call.
sub with_retries {
    my ($what, $code) = @_;
    my $delay = $opt{retry_delay};
    my $msg;
    foreach my $try (0 .. $opt{retries}) {
        (my $r, $msg) = $code->();
        # LWP turns a die inside a request into a failed response
        die("Interrupted by signal $interrupted\n") if $interrupted;
        return $r if defined $r;
        last if ($try == $opt{retries});
        last if (defined $deadline && time() + $delay > $deadline);
        logv("retry ($what): $msg");
        $n{retries}++;
        Time::HiRes::sleep($delay);
        $delay *= 2;
    }
    return (undef, $msg);
}

sub head_once {
    my ($url) = @_;
    my $res = $ua->request(HTTP::Request->new('HEAD', $url));
    return (undef, "HEAD $url: ".$res->status_line) unless $res->is_success;
    return { len => clean_length(scalar($res->header('Content-Length'))),
             etag => clean_etag(scalar($res->header('ETag'))) };
}

# Download to a temporary name beside $dest, verify, set the modification
# time from Last-Modified, sync, rename.  need_lm: fail when the server
# sends no Last-Modified (the modification time of a structure file is
# part of what ParsePDB takes as its identity).
sub fetch_once {
    my ($url, $dest, %o) = @_;
    my $tmp = dirname($dest)."/.".basename($dest).".tmp.$$";
    $TMP{$tmp} = 1;
    my $bad = sub {
        unlink($tmp);
        delete $TMP{$tmp};
        return (undef, "GET $url: $_[0]");
    };
    my $res = $ua->get($url, ':content_file' => $tmp);
    return $bad->($res->status_line) unless $res->is_success;
    return $bad->("transfer aborted: ".($res->header('X-Died') || 'client aborted'))
        if ($res->header('X-Died') || $res->header('Client-Aborted'));
    my $size = -s $tmp;
    return $bad->("empty body") unless $size;
    my $len = clean_length(scalar($res->header('Content-Length')));
    return $bad->("got $size bytes, Content-Length is $len") if (defined $len && $len != $size);
    my $etag = clean_etag(scalar($res->header('ETag')));
    my $md5 = md5_file($tmp);
    my $want = etag_md5($etag);
    return $bad->("MD5 $md5 differs from ETag $want (if the server's ETags are no longer ".
                  "MD5 sums, use --etag-not-md5)") if (defined $want && $want ne $md5);
    # A multipart ETag is checked when it can be.  With one part the
    # part size does not matter, so a difference is an error; with more
    # parts it may only mean a part size that is not known here, and the
    # download is accepted on its length and the gzip test, and counted.
    my $unchecked = 0;
    if (! defined $want && ! $opt{etag_not_md5} && etag_multipart($etag) &&
        ! multipart_matches($tmp, $etag)) {
        return $bad->("the file does not give the one-part ETag $etag (if the server's ETags ".
                      "are no longer made from MD5 sums, use --etag-not-md5)") if ($etag =~ /-1$/);
        $unchecked = 1;
    }
    my $sha = content_sha256($tmp, $dest);
    return $bad->("gzip test failed") unless defined $sha;
    my $lm = $res->header('Last-Modified');
    my $t = defined $lm ? str2time($lm) : undef;
    return $bad->("no usable Last-Modified header") if ($o{need_lm} && ! defined $t);
    if (defined $t) {
        utime($t, $t, $tmp) or return $bad->("cannot set the modification time: $!");
    }
    fsync_path($tmp);
    rename($tmp, $dest) or return $bad->("cannot rename to $dest: $!");
    delete $TMP{$tmp};
    $n{mp_unchecked} += $unchecked;
    return { size => $size, md5 => $md5, sha256 => $sha, etag => $etag,
             mtime => (stat($dest))[9] };
}

# HEAD results of this and of an interrupted run, in the staging
# directory: a rerun does not ask again for a file whose listed date is
# the same.  The line written when a file was fetched also has its MD5
# (fifth column, "-" for the answer to a HEAD request); that is what a
# rerun checks a staged file against when the ETag is not an MD5.
sub load_head_cache {
    my $f = "$stg/HEADS.tsv";
    if (open(my $fh, '<', $f)) {
        while (my $line = <$fh>) {
            next unless ($line =~ s/\n$//);
            my @v = split(/\t/, $line, -1);
            next unless (@v == 5 || @v == 4);
            $headCache{$v[0]} = { listed => $v[1], len => ($v[2] eq '-' ? undef : $v[2]),
                                  etag => ($v[3] eq '-' ? undef : $v[3]) };
            $fetched{$v[0]} = { len => $v[2], etag => $headCache{$v[0]}{etag}, md5 => $v[4] }
                if (@v == 5 && $v[4] =~ /^[0-9a-f]{32}$/);
        }
        close($fh);
    }
    open($headFh, '>>', $f) or die("Cannot write $f: $!\n");
    $headFh->autoflush(1);
}

# $md5: of the file, when this is the record of a download
sub remember_head {
    my ($path, $listed, $head, $md5) = @_;
    $headCache{$path} = { listed => $listed, len => $head->{len}, etag => $head->{etag} };
    $fetched{$path} = { len => $head->{len}, etag => $head->{etag}, md5 => $md5 } if defined $md5;
    print $headFh join("\t", $path, $listed, $head->{len} // '-', $head->{etag} // '-',
                       $md5 // '-'), "\n";
}

sub get_head {
    my ($path, $url, $listed) = @_;
    my $h = $headCache{$path};
    return $h if ($h && $h->{listed} eq $listed);
    my ($head, $msg) = with_retries("HEAD", sub { head_once($url) });
    return (undef, $msg) unless $head;
    remember_head($path, $listed, $head);
    return $head;
}

# Get a structure file into the staging directory.  A file left there
# by an interrupted run is kept only when its bytes are known to be the
# server's: the length is the server's and
#   - the ETag is an MD5 and is the MD5 of the file; or
#   - the ETag is a multipart one and the file gives the same; or
#   - the ETag is the strong ETag the file was fetched with, and the
#     file still has the MD5 it had then.
# Anything else (no ETag, a weak one, a different one that cannot be
# checked) is downloaded again.
sub obtain {
    my ($path, $url, $dest, $listed, $head) = @_;
    $syncDirs{dirname($dest)} = 1;
    if (-f $dest) {
        my $msg;
        ($head, $msg) = get_head($path, $url, $listed) unless $head;
        return (undef, $msg) unless $head;
        my $size = -s $dest;
        if (defined $head->{len} && $size == $head->{len}) {
            my $md5 = md5_file($dest);
            my $want = etag_md5($head->{etag});
            my $was = $fetched{$path};
            my $ok = 0;
            if (defined $want) {
                $ok = ($want eq $md5);
            }
            elsif (defined $head->{etag}) {
                $ok = 1 if (etag_multipart($head->{etag}) && multipart_matches($dest, $head->{etag}));
                $ok = 1 if ($was && defined $was->{etag} && $was->{etag} eq $head->{etag} &&
                            $was->{md5} eq $md5);
            }
            my $sha = $ok ? content_sha256($dest) : undef;
            if (defined $sha) {
                $n{resumed}++;
                return { size => $size, md5 => $md5, sha256 => $sha, etag => $head->{etag},
                         mtime => (stat($dest))[9] };
            }
        }
        unlink($dest) or return (undef, "cannot remove stale $dest: $!");
    }
    make_dir(dirname($dest));
    my ($got, $msg) = with_retries("GET", sub { fetch_once($url, $dest, need_lm => 1) });
    remember_head($path, $listed, { len => $got->{size}, etag => $got->{etag} }, $got->{md5}) if $got;
    return ($got, $msg);
}

# ---------------------------------------------------------------------
# one structure file

sub fail {
    my ($what, $msg) = @_;
    push(@failed, [ $what, $msg ]);
    print STDERR "FAILED $what: $msg\n";
    return undef;
}

# The file an earlier snapshot holds for this path.  $r: its manifest
# row if there is one; otherwise the link in the previous snapshot is
# read.  Returns (file, state, message); state is
#   ok       usable
#   missing  the link or the file does not exist
#   empty    the file has no content
#   size     its size is not the one in the manifest
#   error    anything else (never treated as a missing file)
sub old_file {
    my ($path, $prev, $r) = @_;
    my $storage;
    if ($r) {
        $storage = $r->{storage};
    }
    else {
        my $link = "$data/$prev/snapshot/$path";
        my $l = readlink($link);
        if (! defined $l) {
            return (undef, 'missing') if $!{ENOENT};
            return (undef, 'error', "cannot read the link $link: $!");
        }
        return (undef, 'error', "unexpected link $link -> $l")
            unless ($l =~ m{^(?:\.\./){4}($STORAGE_RE)$});
        $storage = $1;
    }
    my @st = stat("$data/$storage");
    if (! @st) {
        return (undef, 'missing') if $!{ENOENT};
        return (undef, 'error', "cannot stat $data/$storage: $!");
    }
    return (undef, 'error', "$data/$storage is not a plain file") unless (-f _);
    # rec_*: the checksums recorded in the manifest.  They say what the
    # file should be, not what is on disk; see disk_md5.
    my $old = { storage => $storage, size => $st[7], mtime => $st[9], row => $r,
                rec_md5 => ($r && $r->{md5} ne '-' ? $r->{md5} : undef),
                rec_sha256 => ($r && $r->{sha256} ne '-' ? $r->{sha256} : undef) };
    return ($old, 'empty') if ($st[7] == 0);
    return ($old, 'size') if ($r && $r->{size} ne '-' && $r->{size} != $st[7]);
    return ($old, 'ok');
}

# MD5 of the archived file as it is on disk (the file is read, once)
sub disk_md5 {
    my ($old) = @_;
    $old->{disk_md5} = md5_file("$data/$old->{storage}") unless defined $old->{disk_md5};
    return $old->{disk_md5};
}

# SHA-256 of its content as it is on disk; undef when it fails the gzip test
sub disk_sha256 {
    my ($old) = @_;
    $old->{disk_sha256} = content_sha256("$data/$old->{storage}")
        unless exists $old->{disk_sha256};
    return $old->{disk_sha256};
}

# MD5 to compare with what the server says before anything is
# downloaded: the recorded one, so that an unchanged file is not read;
# the file is read when nothing is recorded.  Whether the file still
# has the recorded checksum is what --verify and --refetch look at.
sub old_md5 {
    my ($old) = @_;
    return $old->{rec_md5} // disk_md5($old);
}

# Does the file on disk have the checksums recorded for it?  1 yes,
# 0 no, undef when nothing is recorded.  The MD5 settles it when there
# is one; $both also decompresses the file and compares the SHA-256.
sub agrees_with_record {
    my ($old, $both) = @_;
    return undef unless (defined $old->{rec_md5} || defined $old->{rec_sha256});
    return 0 if (defined $old->{rec_md5} && disk_md5($old) ne $old->{rec_md5});
    if (defined $old->{rec_sha256} && ($both || ! defined $old->{rec_md5})) {
        my $sha = disk_sha256($old);
        return 0 unless (defined $sha && $sha eq $old->{rec_sha256});
    }
    return 1;
}

sub make_links {
    my ($k, $e, $storage) = @_;
    my $name = entry_name($k, $e);
    my $hh = substr($e, 1, 2);
    my $hdir = "$stg/snapshot/$k->{dir}/$hh";
    make_dir($hdir);
    $syncDirs{$hdir} = 1;
    symlink("../../../../$storage", "$hdir/$name")
        or die("Cannot make link $hdir/$name: $!\n");
    my $all = "$stg/snapshot/".all_dir($k)."/$name";
    symlink("../$k->{dir}/$hh/$name", $all) or die("Cannot make link $all: $!\n");
}

# Returns the manifest row, or undef after recording a failure.
#   $state  new, cand (HEAD wanted) or keep
#   $r      row of the previous manifest for this path, if any
#   $rp     set when an earlier snapshot listed this (new) file as removed
sub process_entry {
    my ($k, $e, $state, $listed, $prev, $r, $rp, $reappeared) = @_;
    my $kind = $k->{kind};
    my $path = member_path($k, $e);
    my $url = "$opt{source}/".remote_path($k, $e);
    my $dest = "$stg/files/$path";
    my $failed = sub { $c{$kind}{failed}++; return fail($path, $_[0]); };

    # the copy already in the archive, if any
    my ($old, $ost, $err) = (undef, 'none');
    if ($state ne 'new') {
        ($old, $ost, $err) = old_file($path, $prev, $r);
    }
    elsif ($rp) {
        $r = $rp->{row};
        ($old, $ost, $err) = old_file($path, undef, $r) if ($r->{storage} ne '-');
    }
    return $failed->($err) if ($ost eq 'error');
    # --verify: the file must still have the checksums recorded for it
    if ($opt{verify} && $ost eq 'ok') {
        my $agrees = agrees_with_record($old, 1);
        if (! defined $agrees) {
            $n{unverified}++;
        }
        elsif ($agrees) {
            $n{verified}++;
        }
        else {
            $ost = 'checksum';
        }
    }
    if ($state eq 'new') {
        # a kept copy of a file that comes back is used only if it is good
        ($old, $ost) = (undef, 'none') if ($ost ne 'ok');
        $state = 'cand' if $old;
    }

    # A missing copy is replaced from the server, within limits.  A copy
    # that is empty, has the wrong size or (--verify) the wrong checksum
    # fails the file unless --refetch names it.  --refetch always
    # downloads, and the download is compared with the bytes on disk.
    my ($repair, $forced);
    if ($ost eq 'missing') {
        $repair = 'missing';
        $n{repaired}++;
        return $failed->("the archived file is missing, and more than $opt{max_repaired} are")
            if ($n{repaired} > $opt{max_repaired});
    }
    elsif ($ost eq 'empty' || $ost eq 'size' || $ost eq 'checksum') {
        my $what = ($ost eq 'empty') ? "is empty"
            : ($ost eq 'size') ? "is $old->{size} bytes, its manifest says $r->{size}"
            : "does not have the checksums recorded in its manifest";
        return $failed->("the archived file $old->{storage} $what.  To replace it from the ".
                         "server after comparing the content, run again with --refetch $path")
            unless $refetch{$path};
        if ($ost eq 'empty') {
            $repair = 'empty';
        }
        else {
            $forced = $ost;
        }
    }
    elsif ($old && $refetch{$path}) {
        $forced = 'asked';
    }
    $refetch{$path} = 2 if ($forced || ($repair && $repair eq 'empty'));
    # nothing to compare a download with
    $old = undef if $repair;

    # is a download needed?
    my ($head, $why);
    if (! $old) {
        $why = $repair ? 'repair' : 'new';
    }
    elsif ($forced) {
        $why = 'refetch';
    }
    elsif ($state eq 'cand') {
        ($head, my $msg) = get_head($path, $url, $listed);
        return $failed->($msg) unless $head;
        $c{$kind}{checked}++;
        my $want = etag_md5($head->{etag});
        my $sameLen = (defined $head->{len} && $head->{len} == $old->{size});
        if (defined $want) {
            $why = 'changed' unless ($sameLen && $want eq old_md5($old));
        }
        elsif ($sameLen && etag_multipart($head->{etag}) &&
               (($r && $r->{etag} eq $head->{etag}) ||
                multipart_matches("$data/$old->{storage}", $head->{etag}))) {
            # multipart ETag: the one recorded for this file, or the
            # one the local file gives
        }
        else {
            # no usable identity from the server: compare the content itself
            $why = 'unknown';
        }
    }

    my $got;
    if ($why) {
        ($got, my $msg) = obtain($path, $url, $dest, $listed, $head);
        return $failed->($msg) unless $got;
        $head = { etag => $got->{etag} };
        if ($old) {
            # Compared with the archived file as it is on disk, never
            # with what the manifest says about it.  The same bytes, or
            # the same content compressed differently: keep that file.
            my $same;
            if ($got->{md5} eq disk_md5($old)) {
                $same = 'identical';
                $old->{disk_sha256} = $got->{sha256};
            }
            elsif (defined disk_sha256($old) && $got->{sha256} eq $old->{disk_sha256}) {
                $same = 'recompressed';
            }
            # does the archived file have the size and checksum recorded for it?
            my $agrees = ($ost eq 'size' || $ost eq 'checksum') ? 0 : agrees_with_record($old);
            if ($same) {
                unlink($dest) or die("Cannot remove $dest: $!\n");
                $c{$kind}{$same}++;
                $n{notkept}++;
                $got = undef;
                my $wrong = ($ost eq 'size') ? "size" : (defined $agrees && ! $agrees) ? "checksum" : undef;
                push(@notes, "$path: archived file agrees with the server, kept".
                     ($wrong ? " (manifest $wrong was wrong; the new manifest has the file's own)" : ""))
                    if ($forced || $wrong);
            }
            elsif ($forced && ! $agrees) {
                # damaged, or nothing recorded to tell: the download replaces it
                $repair = ($forced ne 'asked') ? $forced : defined $agrees ? 'checksum' : 'asked';
            }
            elsif (defined $agrees && ! $agrees && ! $rp) {
                # found damaged in passing; the download waits in the
                # staging directory for the run with --refetch
                return $failed->("the archived file $old->{storage} does not have the checksum ".
                                 "recorded in its manifest, and is not what the server has.  To ".
                                 "replace it from the server, run again with --refetch $path");
            }
            else {
                # the archived file is as recorded and the server has
                # something else
                $why = 'changed' if $forced;
            }
        }
    }

    my %row = (path => $path, entry => $e, kind => $kind,
               listed_date => ($listed eq '' ? '-' : $listed));
    if ($got) {
        if ($repair) {
            # the same content as recorded for the lost file gets its
            # modification time back
            my $hasRec = ($r && ($r->{md5} ne '-' || $r->{sha256} ne '-'));
            my $known = ($r && (($r->{md5} ne '-' && $r->{md5} eq $got->{md5}) ||
                                ($r->{sha256} ne '-' && $r->{sha256} eq $got->{sha256})));
            my $t = $known ? iso_epoch($r->{last_modified}) : undef;
            my $was = { missing => "was missing", empty => "was empty",
                        size => "had the wrong size",
                        checksum => "did not have its recorded checksum",
                        asked => "differs from the server's (--refetch)" }->{$repair};
            if (defined $t) {
                if (! utime($t, $t, $dest)) {
                    unlink($dest);
                    return $failed->("cannot set the modification time of $dest: $!");
                }
                $got->{mtime} = $t;
                push(@notes, "$path: archived file $was, fetched again; content and ".
                     "modification time as recorded");
            }
            else {
                push(@notes, "$path: archived file $was, fetched again; ".
                     ($hasRec && ! $known ? "the download is not the recorded content either"
                      : $known ? "no recorded modification time"
                      : "no record to compare with").
                     ", so ParsePDB may take it as a new release");
            }
            $c{$kind}{repaired}++;
        }
        my $storage = "$date/files/$path";
        make_links($k, $e, $storage);
        # a file that comes back is new to this snapshot, whatever its content
        $c{$kind}{ ($why eq 'new' || $rp) ? 'new' : $repair ? 'refetched' : 'changed' }++;
        logv("fetched ($why) $path");
        @row{qw(storage size md5 sha256 last_modified etag status source)} =
            ($storage, $got->{size}, $got->{md5}, $got->{sha256}, iso_time($got->{mtime}),
             $got->{etag} // '-', 'fetched', $date);
    }
    else {
        make_links($k, $e, $old->{storage});
        $c{$kind}{reused}++;
        logv("reused $path -> $old->{storage}");
        # what was read from the file itself in this run, else the record
        my $md5 = $old->{disk_md5} // $old->{rec_md5} // '-';
        my $sha = $old->{disk_sha256} // $old->{rec_sha256} // '-';
        if ($opt{hash_reused}) {
            $md5 = disk_md5($old) if ($md5 eq '-');
            $sha = disk_sha256($old) if ($sha eq '-');
            return $failed->("$old->{storage} fails the gzip test") unless defined $sha;
        }
        my $etag = ($head && defined $head->{etag}) ? $head->{etag} : ($r ? $r->{etag} : '-');
        @row{qw(storage size md5 sha256 last_modified etag status source)} =
            ($old->{storage}, $old->{size}, $md5, $sha, iso_time($old->{mtime}),
             $etag, 'reused', substr($old->{storage}, 0, 6));
    }
    if ($rp) {
        push(@$reappeared, [ $path, $rp->{snap}, ($got ? 'fetched' : "reused $old->{storage}") ]);
        push(@notes, "$path: removed in $rp->{snap}, listed again; ".
             ($got ? "fetched (no identical archived copy)"
                   : "identical to $old->{storage}, which is linked; ParsePDB will not see ".
                     "it in files/ of this snapshot"));
    }
    return \%row;
}

# ---------------------------------------------------------------------
# sequences and index files

sub aux_row {
    my ($rel, $kind, $got) = @_;
    my $f = "$stg/files/$rel";
    fsync_path($f);
    my @st = stat($f);
    return { path => $rel, entry => '-', kind => $kind, storage => "$date/files/$rel",
             size => $st[7], md5 => $got->{md5} // md5_file($f),
             sha256 => $got->{sha256} // content_sha256($f),
             last_modified => iso_time($st[9]), etag => $got->{etag} // '-',
             listed_date => '-', status => 'fetched', source => $date };
}

sub get_seqres {
    my $rel = "seq/pdbseqres.$date.fa";
    my $out = "$stg/files/$rel";
    my $url = "$opt{source}/$opt{seqres_path}";
    my $gz = ($url =~ /\.gz$/);
    my $dl = $gz ? "$out.gz" : $out;
    my ($got, $msg) = with_retries("GET", sub { fetch_once($url, $dl) });
    return fail($rel, $msg) unless $got;
    if ($gz) {
        # stored uncompressed, as before
        my $tmp = "$out.tmp.$$";
        $TMP{$tmp} = 1;
        my $in = open_cmd($opt{gzip}, '-dc', '--', $dl);
        open(my $o, '>', $tmp) or die("Cannot write $tmp: $!\n");
        binmode($o);
        my $buf;
        while (read($in, $buf, 1 << 20)) {
            print $o $buf or die("Cannot write $tmp: $!\n");
        }
        close($o) or die("Cannot write $tmp: $!\n");
        close($in) or return fail($rel, "cannot decompress $dl");
        utime($got->{mtime}, $got->{mtime}, $tmp)
            or return fail($rel, "cannot set the modification time: $!");
        rename($tmp, $out) or die("Cannot rename $tmp to $out: $!\n");
        delete $TMP{$tmp};
        unlink($dl);
        $got = { sha256 => $got->{sha256}, etag => $got->{etag} };
    }
    return fail($rel, "file is empty") unless (-s $out);
    open(my $fh, '<', $out) or die("Cannot read $out: $!\n");
    my $first = <$fh>;
    close($fh);
    return fail($rel, "does not start with a FASTA header") unless (defined $first && $first =~ /^>/);
    symlink("pdbseqres.$date.fa", "$stg/files/seq/pdbseqres.fa")
        or die("Cannot make link pdbseqres.fa: $!\n");
    return aux_row($rel, 'seq', $got);
}

sub get_onhold {
    my $rel = "seq/onhold.$date.fa";
    my $out = "$stg/files/$rel";
    my ($got, $msg) = with_retries("GET", sub {
        my ($g, $m) = fetch_once($opt{onhold_url}, $out);
        return (undef, $m) unless $g;
        # must be FASTA: first non-blank line a header, at least one sequence
        open(my $fh, '<', $out) or die("Cannot read $out: $!\n");
        my ($first, $headers, $seqs) = (undef, 0, 0);
        while (my $l = <$fh>) {
            next if ($l =~ /^\s*$/);
            $first = $l unless defined $first;
            if ($l =~ /^>/) { $headers++; }
            elsif ($l =~ /^[A-Za-z*-]+\s*$/) { $seqs++; }
        }
        close($fh);
        if (! defined $first || $first !~ /^>/ || ! $headers || ! $seqs) {
            unlink($out);
            return (undef, "on-hold file is not FASTA (URL changed?)");
        }
        return $g;
    });
    return fail($rel, $msg) unless $got;
    symlink("onhold.$date.fa", "$stg/files/seq/onhold.fa")
        or die("Cannot make link onhold.fa: $!\n");
    return aux_row($rel, 'seq', $got);
}

# protein sequences only, sorted by header
sub make_pdbprot {
    my $rel = "seq/pdbprot.$date.fa";
    my $out = "$stg/files/$rel";
    my %seq;
    {
        local $/ = "\n>";
        open(my $in, '<', "$stg/files/seq/pdbseqres.$date.fa")
            or die("Cannot read pdbseqres.$date.fa: $!\n");
        while (<$in>) {
            my ($header, $seq) = /^>?([^\n]*)\n(.*)$/s;
            next unless defined $header;
            if ($header =~ /\smol\:protein/) { # both protein and protein-het
                $seq =~ s/[\n>]//sg;
                $seq{$header} = $seq;
            }
        }
        close($in);
    }
    return fail($rel, "no protein sequences in pdbseqres.$date.fa") unless %seq;
    open(my $o, '>', $out) or die("Cannot write $out: $!\n");
    foreach my $header (sort(keys(%seq))) {
        print $o ">$header\n".$seq{$header}."\n" or die("Cannot write $out: $!\n");
    }
    close($o) or die("Cannot write $out: $!\n");
    symlink("pdbprot.$date.fa", "$stg/files/seq/pdbprot.fa")
        or die("Cannot make link pdbprot.fa: $!\n");
    return aux_row($rel, 'seq', {});
}

sub get_index {
    my $base = "$opt{source}/$opt{index_path}";
    $base .= "/" unless ($base =~ m{/$});
    my ($html, $msg) = with_retries("GET", sub {
        my $res = $ua->get($base);
        return (undef, "GET $base: ".$res->status_line) unless $res->is_success;
        return $res->decoded_content;
    });
    return fail("index", $msg) unless defined $html;
    # file names: relative links, or links to this same directory
    my %names;
    foreach my $href ($html =~ /<a\s[^>]*href="([^"]+)"/gi) {
        $href =~ s/^\Q$base\E//;
        next if ($href =~ m{[/?#]} || $href =~ /^\.\.?$/);
        $href =~ s/%([0-9A-Fa-f]{2})/chr(hex($1))/ge;
        next unless ($href =~ /^[A-Za-z0-9][A-Za-z0-9._+-]*$/);
        $names{$href} = 1;
    }
    my @missing = grep { ! $names{$_} } split(/,/, $opt{index_required});
    return fail("index", "listing $base lacks: ".join(", ", @missing)) if @missing;
    # every file is stored gzipped; one that is listed gzipped is kept as it is
    my @clash = grep { /\.gz$/ && $names{substr($_, 0, -3)} } keys %names;
    return fail("index", "listing $base has both a file and its .gz: ".join(", ", sort @clash))
        if @clash;
    my @rows;
    foreach my $name (sort keys %names) {
        my $f = "$stg/files/index/$name";
        my ($got, $m) = with_retries("GET", sub { fetch_once("$base$name", $f) });
        return fail("index/$name", $m) unless $got;
        if ($name =~ /\.gz$/) {
            push(@rows, aux_row("index/$name", 'index', $got));
            next;
        }
        unlink("$f.gz");
        if (system({ $opt{gzip} } $opt{gzip}, '-9', '--', $f) != 0 || ! -s "$f.gz") {
            return fail("index/$name", "gzip failed");
        }
        push(@rows, aux_row("index/$name.gz", 'index', { sha256 => $got->{sha256}, etag => $got->{etag} }));
    }
    return @rows;
}

# ---------------------------------------------------------------------

sub summary {
    info("");
    info(sprintf("%-8s %8s %8s %8s %8s %8s %8s", "kind", "fetched", "reused",
                 "checked", "removed", "failed", "repaired"));
    my %t;
    foreach my $k (@KINDS) {
        my $x = $c{$k->{kind}} || {};
        my $fetched = ($x->{new} || 0) + ($x->{changed} || 0) + ($x->{refetched} || 0);
        info(sprintf("%-8s %8d %8d %8d %8d %8d %8d", $k->{kind}, $fetched,
                     $x->{reused} || 0, $x->{checked} || 0, $x->{removed} || 0,
                     $x->{failed} || 0, $x->{repaired} || 0));
        $t{fetched} += $fetched;
        $t{$_} += ($x->{$_} || 0) foreach (qw(changed reused removed failed identical recompressed repaired));
    }
    info("total: fetched $t{fetched} (of which changed content $t{changed}), ".
         "reused $t{reused}, removed $t{removed}, failed ".scalar(@failed));
    info("downloaded but identical to the archived file, not kept: $t{identical}") if $t{identical};
    info("downloaded, same content compressed differently, not kept: $t{recompressed}") if $t{recompressed};
    info("damaged or missing archived files fetched again: $t{repaired}") if $t{repaired};
    info("verified against the previous manifest (--verify): $n{verified} archived files; ".
         "$n{unverified} have no recorded checksum and were not verified") if $opt{verify};
    info("downloads with a multipart ETag that could not be checked against it (part size ".
         "not known; accepted on length and gzip test): $n{mp_unchecked}") if $n{mp_unchecked};
    info("kept from an interrupted run: $n{resumed}") if $n{resumed};
    info("retries: $n{retries}") if $n{retries};
    info("note: ETags were not used as MD5 sums (--etag-not-md5)") if $opt{etag_not_md5};
    info("note: --refetch $_ did nothing: no archived file of that name was to be reused ".
         "(the file is new, not listed, or its archived copy is missing)")
        foreach (grep { $refetch{$_} == 1 } sort keys %refetch);
    info("note: $_") foreach (@notes);
    info("failed: $_->[0]: $_->[1]") foreach (@failed);
}

# ---------------------------------------------------------------------
# --manifest-only: write or complete the manifest of an existing snapshot.
# Reads every file that has no checksum yet; with --verify, every file,
# and the checksums the manifest already has must be those of the file.
# A file that is read is always compared with what is recorded for it.
# Can be stopped and run again: finished rows are appended to
# <manifest>.partial as they are made.

# A snapshot whose directory is not writable by its owner is frozen (it is
# the freeze of a SCOP/SCOPe release) and is never changed, whoever runs
# this; -w alone would let root through.
sub frozen {
    my ($d) = @_;
    my @st = stat($d);
    return (! @st) || (! ($st[2] & 0200)) || (! -w $d);
}

sub manifest_only {
    my ($snap, $reference) = @_;
    die("'$snap' is not a snapshot name (YYMMDD)\n") unless defined snapshot_day($snap);
    my $sdir = "$data/$snap";
    die("Snapshot $sdir does not exist\n") unless (-d "$sdir/files");
    my $own = "$sdir/files/MANIFEST.tsv";
    my $final = $opt{manifest_out} // $own;
    my $partial = "$final.partial";
    if (frozen($sdir)) {
        my $outDir = defined $opt{manifest_out} ? realpath(dirname($opt{manifest_out})) : undef;
        my $real = realpath($sdir);
        die("Snapshot $snap is frozen (read-only) and is never written to; "
            ."give --manifest-out FILE outside it\n")
            if (! defined $outDir || $outDir eq $real || index("$outDir/", "$real/") == 0);
    }
    die("--reference cannot be used with --verify: every file is read\n")
        if ($opt{verify} && @$reference);

    # rows wanted: {path} = line
    my ($rows, $prevName);
    if (-f $own) {
        $rows = read_manifest($own);
        open(my $fh, '<', $own) or die("Cannot read $own: $!\n");
        while (<$fh>) {
            last unless /^#/;
            $prevName = $1 if /^#previous\t(\S+)/;
        }
        close($fh);
        info("completing $own");
    }
    else {
        die("Snapshot $snap is pruned; run makelinks.pl first\n") unless (-d "$sdir/snapshot");
        info("building the manifest of $snap from its link tree");
        my $broken = 0;
        $rows = {};
        foreach my $k (@KINDS) {
            my $d = "$sdir/snapshot/".all_dir($k);
            opendir(my $dh, $d) or die("Cannot read $d: $!\n");
            my @names = grep { $_ ne '.' && $_ ne '..' } readdir($dh);
            closedir($dh);
            foreach my $name (@names) {
                my $e = name_entry($k, $name);
                next unless defined $e;
                my $path = member_path($k, $e);
                my $l = readlink("$sdir/snapshot/$path");
                if (! defined $l || $l !~ m{^(?:\.\./){4}($STORAGE_RE)$}) {
                    print STDERR "broken link: $snap/snapshot/$path\n";
                    $broken++;
                    next;
                }
                my $storage = $1;
                my $src = substr($storage, 0, 6);
                $rows->{$path} = join("\t", $path, $e, $k->{kind}, $storage, ('-') x 6,
                                      ($src eq $snap ? 'fetched' : 'reused'), $src);
            }
        }
        foreach my $sub ('seq', 'index') {
            my $d = "$sdir/files/$sub";
            next unless (-d $d);
            opendir(my $dh, $d) or die("Cannot read $d: $!\n");
            foreach my $name (grep { $_ !~ /^\./ } readdir($dh)) {
                next if (-l "$d/$name" || ! -f "$d/$name");
                $rows->{"$sub/$name"} = join("\t", "$sub/$name", '-', $sub,
                                             "$snap/files/$sub/$name", ('-') x 6, 'fetched', $snap);
            }
            closedir($dh);
        }
        die("$broken broken link(s) in $snap; no manifest written\n") if $broken;
    }

    # checksums already known, by storage path: [size, md5, sha256]
    my %known;
    open(my $pfh, '+>>', $partial) or die("Cannot open $partial: $!\n");
    flock($pfh, LOCK_EX | LOCK_NB) or die("Another run is writing $partial\n");
    foreach my $f ($partial, @$reference) {
        my $ref = read_manifest($f, ($f eq $partial));
        foreach my $line (values %$ref) {
            my @v = split(/\t/, $line, -1);
            next if ($v[4] eq '-' || $v[5] eq '-' || $v[6] eq '-');
            $known{$v[3]} = join("\t", @v[4, 5, 6]);
        }
    }
    seek($pfh, 0, 2);
    $pfh->autoflush(1);

    my ($hashed, $taken, $bad, $done, $filled, $compared, $differ) = (0, 0, 0, 0, 0, 0, 0);
    my $total = scalar(keys %$rows);
    my @out;
    foreach my $path (sort keys %$rows) {
        $done++;
        info("$done of $total rows; $hashed files read") if ($done % 10000 == 0);
        my $r = row_hash($rows->{$path});
        my $f = "$data/$r->{storage}";
        my @st = stat($f);
        if (! @st || ! -f _) {
            print STDERR "missing: $r->{storage}\n";
            $bad++;
            next;
        }
        if ($r->{size} ne '-' && $r->{size} != $st[7]) {
            print STDERR "size differs from manifest: $r->{storage}\n";
            $bad++;
            next;
        }
        my $missing = ($r->{md5} eq '-' || $r->{sha256} eq '-');
        $filled++ if ($missing || $r->{size} eq '-' || $r->{last_modified} eq '-');
        $r->{size} = $st[7];
        $r->{last_modified} = iso_time($st[9]) if ($r->{last_modified} eq '-');
        if ($missing || $opt{verify}) {
            # Checksums from the work file (this file was read by an
            # interrupted pass) or from a reference manifest are taken
            # only when they do not contradict what is recorded here.
            my @k = split(/\t/, $known{$r->{storage}} // '');
            if (@k && $k[0] == $st[7] && ($r->{md5} eq '-' || $r->{md5} eq $k[1]) &&
                ($r->{sha256} eq '-' || $r->{sha256} eq $k[2])) {
                @$r{qw(md5 sha256)} = @k[1, 2];
                $taken++;
            }
            else {
                my $sha = content_sha256($f);
                if (! defined $sha) {
                    print STDERR "fails gzip test: $r->{storage}\n";
                    $bad++;
                    next;
                }
                my $md5 = md5_file($f);
                $hashed++;
                $compared++ if ($r->{md5} ne '-' || $r->{sha256} ne '-');
                my @off = grep { $r->{$_->[0]} ne '-' && $r->{$_->[0]} ne $_->[1] }
                    ([ 'md5', $md5 ], [ 'sha256', $sha ]);
                if (@off) {
                    print STDERR "checksum differs from manifest (".join(", ", map { $_->[0] } @off).
                        "): $r->{storage}\n";
                    $bad++;
                    $differ++;
                    next;
                }
                $r->{md5} = $md5;
                $r->{sha256} = $sha;
                print $pfh join("\t", map { $r->{$_} } @COLS), "\n";
            }
        }
        push(@out, join("\t", map { $r->{$_} } @COLS));
    }
    info("$total rows: $hashed files read, $taken checksums taken from other manifests, $bad bad");
    info("verify: $compared files read and compared with their recorded checksums, $differ differ",
         ($taken ? "; $taken taken as read by an interrupted pass ($partial)" : ""))
        if $opt{verify};
    if ($bad) {
        info("No manifest written; finished rows are kept in $partial");
        return 1;
    }
    if ($opt{verify} && ! $filled && $final eq $own && -f $own) {
        close($pfh);
        unlink($partial);
        info("$own was complete and is unchanged");
        return 0;
    }
    # in the order snapshot.pl writes: by kind, then path
    my %order;
    my $i = 0;
    $order{$_->{kind}} = $i++ foreach (@KINDS);
    $order{seq} = $i++;
    $order{index} = $i++;
    my @kindOf = map { $order{(split(/\t/, $_, 4))[2]} } @out;
    @out = @out[ sort { $kindOf[$a] <=> $kindOf[$b] || $out[$a] cmp $out[$b] } 0 .. $#out ];
    write_manifest($final, $snap, $prevName, \@out);
    close($pfh);
    unlink($partial);
    info("wrote $final");
    return 0;
}
