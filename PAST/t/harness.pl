#!/usr/bin/perl
#
#   Copyright (C) 2026 John-Marc Chandonia
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
#   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
#   02110-1301, USA.
#
#   This file is part of the PAST (PDB Archival Snapshot Toolkit) package.
#
#
# Test harness for PAST snapshot.pl, prune.pl and makelinks.pl.  Builds a
# small fake archive and a fake wwPDB site served by t/mockserver.py, and
# runs the scripts against them.  Nothing outside the work directory is
# touched.
#
# usage: perl t/harness.pl          (or t/run)
# environment:
#   PAST_TEST_TMP    where to make the work directory (default: system tmp)
#   PAST_TEST_KEEP   keep the work directory
#   PAST_TEST_SHOW   print the output of every run
#
# snapshot.pl runs with its default guards and limits unless a case says
# otherwise; only the retry delay and the number of retries are reduced.

use strict;
use warnings;
use Test::More;
use File::Temp qw(tempdir);
use File::Path qw(make_path remove_tree);
use File::Basename qw(dirname);
use File::Find;
use Cwd qw(abs_path);
use IO::Compress::Gzip qw(gzip $GzipError);
use Digest::MD5 qw(md5_hex);
use Digest::SHA qw(sha256_hex);
use Fcntl qw(:flock);
use JSON::PP;
use POSIX ();
use Time::HiRes qw(sleep);
use Time::Local qw(timegm);

my $here = dirname(abs_path($0));
my $snapshot = $ENV{PAST_SNAPSHOT} // "$here/../snapshot.pl";
my $prune = $ENV{PAST_PRUNE} // "$here/../prune.pl";
my $makelinks = $ENV{PAST_MAKELINKS} // "$here/../makelinks.pl";

my $work = tempdir("past-test.XXXXXX", DIR => ($ENV{PAST_TEST_TMP} // File::Spec->tmpdir),
                   CLEANUP => ($ENV{PAST_TEST_KEEP} ? 0 : 1));
my ($root, $doc, $ctl, $tmp) = ("$work/root", "$work/docroot", "$work/control", "$work/tmp");
make_path($doc, $ctl, $tmp);
my $pub = "$doc/pub/pdb";

# ---------------------------------------------------------------- site

my %site;       # id => { obsolete, date, files => { key => text }, path => { key => p } }
my %hold;       # options that damage the holdings
my %lm;         # remote path => mtime to set
my %gzLevel;    # remote path => compression level (default 6)
my @indexNames = qw(author.idx cmpd_res.idx compound.idx crystal.idx entries.idx
                    on_hold.list pending_waiting.list resolu.idx source.idx
                    status_query.csv status_query.seq);
my %indexHidden;
my %indexExtra; # name => bytes
my $onhold = ">1XYZ_1 mol:protein length:5  ON HOLD\nMKVLA\n>1XYW_1 mol:protein length:4  ON HOLD\nGGSA\n";

# snapshot dates are a week apart; listed dates fall between them
my $t0 = timegm(0, 0, 0, 1, 0, 2030);
my $week = 0;
my $D;          # date of the run being made
my $L;          # last published snapshot
sub snap_date { return POSIX::strftime("%y%m%d", gmtime($t0 + 7 * 86400 * $_[0])); }
sub next_date { $week++; $D = snap_date($week); return $D; }
sub listed_now {
    return POSIX::strftime("%Y-%m-%dT10:00:00+00:00", gmtime($t0 + 7 * 86400 * $week + 4 * 86400));
}

sub rpath {
    my ($id, $key, $obs) = @_;
    my $hh = substr($id, 1, 2);
    my $d = $obs ? "obsolete" : "divided";
    return "data/structures/$d/pdb/$hh/pdb$id.ent.gz"   if ($key eq 'pdb');
    return "data/structures/$d/mmCIF/$hh/$id.cif.gz"    if ($key eq 'mmcif');
    return "data/structures/$d/XML/$hh/$id.xml.gz"      if ($key eq 'pdbml');
    return "compatible/pdb_bundle/$hh/$id/$id-pdb-bundle.tar.gz" if ($key eq 'pdb_bundle');
    die;
}

sub gz {
    my ($text, $level) = @_;
    my $out;
    gzip(\$text => \$out, Time => 0, Minimal => 1, Level => ($level // 6)) or die $GzipError;
    return $out;
}

sub put {
    my ($file, $bytes) = @_;
    if (-f $file) {
        open(my $fh, '<', $file) or die;
        binmode($fh);
        local $/;
        my $old = <$fh>;
        return if ($old eq $bytes);
    }
    make_path(dirname($file));
    open(my $fh, '>', $file) or die "$file: $!";
    binmode($fh);
    print $fh $bytes;
    close($fh);
}

sub add_entry {
    my ($id, %o) = @_;
    my @keys = $o{keys} ? @{$o{keys}} : qw(pdb mmcif pdbml);
    $site{$id} = { obsolete => $o{obsolete} || 0,
                   date => (exists $o{date} ? $o{date} : listed_now()),
                   files => { map { $_ => "$_ file of $id, version 1\n".
                                    join("", map { "filler line $_ of a structure file\n" } 1 .. 40) } @keys } };
}

sub bump { $site{$_}{date} = listed_now() foreach (@_); }

# number of structure files a snapshot of the site must have
sub nstruct {
    my $n = 0;
    foreach my $s (values %site) {
        $n += scalar(grep { ! $s->{path}{$_} } keys %{$s->{files}});
    }
    return $n;
}

sub write_site {
    my (%cur, %obs, %curDates, %obsDates, %wanted);
    foreach my $id (sort keys %site) {
        my $s = $site{$id};
        my %ct;
        foreach my $key (keys %{$s->{files}}) {
            my $p = $s->{path}{$key} // rpath($id, $key, $s->{obsolete});
            $wanted{$p} = 1;
            put("$pub/$p", gz($s->{files}{$key}, $gzLevel{$p}));
            utime($lm{$p}, $lm{$p}, "$pub/$p") if $lm{$p};
            my $hkey = ($key eq 'pdb' && $hold{rename_pdb_key}) ? $hold{rename_pdb_key} : $key;
            $ct{$hkey} = [ "/pdb/$p" ];
        }
        if ($s->{obsolete}) {
            $obs{uc $id} = { status_code => "OBS", content_type => \%ct };
            $obsDates{uc $id} = $s->{date} if defined $s->{date};
        }
        else {
            $cur{uc $id} = \%ct;
            $curDates{uc $id} = $s->{date};
        }
    }
    # files no longer listed disappear from the site
    foreach my $top ("$pub/data", "$pub/compatible") {
        next unless (-d $top);
        find(sub { return unless -f $_;
                   (my $rel = $File::Find::name) =~ s{^\Q$pub\E/}{};
                   unlink($_) unless $wanted{$rel}; }, $top);
    }
    my $json = JSON::PP->new->canonical;
    put("$pub/holdings/current_file_holdings.json.gz", gz($json->encode(\%cur)));
    put("$pub/holdings/all_removed_entries.json.gz", gz($json->encode(\%obs)));
    put("$pub/holdings/released_structures_last_modified_dates.json.gz", gz($json->encode(\%curDates)));
    put("$pub/holdings/obsolete_structures_last_modified_dates.json.gz", gz($json->encode(\%obsDates)));

    put("$pub/derived_data/pdb_seqres.txt.gz", gz(
        ">1aaa_A mol:protein length:6  PROTEIN ONE\nMKVLAA\n".
        ">1aaa_B mol:na length:4  DNA\nACGT\n".
        ">1aab_A mol:protein length:5  PROTEIN TWO\nGGSAW\n"));
    # the listing as files.wwpdb.org sends it: relative links, a parent
    # link and a link to another site
    my $html = '<!DOCTYPE html><html><body><table><tbody><tr><td colspan="3">'.
        '<a href="../">Parent Directory</a></td></tr>';
    foreach my $n (@indexNames, sort keys %indexExtra) {
        put("$pub/derived_data/index/$n", $indexExtra{$n} // "contents of $n\n");
        next if $indexHidden{$n};
        $html .= '<tr><td><i class="bi me-2"></i><a href="'.$n.'">'.$n.'</a></td><td>2029-12-14 12:27:21</td></tr>';
    }
    $html .= '</tbody></table><p>See the <a href="https://www.wwpdb.org/ftp/pdb-ftp-sites" target="_blank">wwPDB website</a>.</p></body></html>';
    put("$pub/derived_data/index/index.html", $html);
    put("$doc/onhold", $onhold);
}

my %faults;
sub set_faults {
    %faults = @_;
    # the on-hold service sends no Last-Modified and a weak ETag
    $faults{"/onhold"} //= { no_last_modified => JSON::PP::true, etag_raw => 'W/"33980-Hqqh698"' };
    open(my $fh, '>', "$ctl/faults.json.tmp") or die;
    print $fh JSON::PP->new->encode(\%faults);
    close($fh);
    rename("$ctl/faults.json.tmp", "$ctl/faults.json") or die;
}
sub url_of { return "/pub/pdb/".rpath(@_); }

# ---------------------------------------------------------------- helpers

my $base;
my $lastOut = "";

sub shq { my ($s) = @_; $s =~ s/'/'\\''/g; return "'$s'"; }

sub snap_cmd {
    return ("timeout", "120", $^X, $snapshot, "--root", $root, "--source", "$base/pub/pdb",
            "--tmpdir", $tmp, "--onhold-url", "$base/onhold", "--retries", 2,
            "--retry-delay", "0.05", @_);
}

# runs snapshot.pl without a terminal on standard input; returns exit status
sub snap {
    my $cmd = join(" ", map { shq($_) } snap_cmd(@_));
    $lastOut = `$cmd </dev/null 2>&1`;
    my $rc = $? >> 8;
    print STDERR "\n--- snapshot.pl @_ (exit $rc)\n$lastOut---\n" if $ENV{PAST_TEST_SHOW};
    return $rc;
}

# a run that must publish $D
sub publish {
    my ($label, @args) = @_;
    my $ok = is(snap("--date", $D, @args), 0, $label);
    $ok = is(latest(), $D, "latest is $D") && $ok;
    diag($lastOut) unless $ok;
    $L = $D if $ok;
    return $ok;
}

# a run that must fail and leave the archive as it was
sub refused {
    my ($label, $re, @args) = @_;
    my $before = ls_data();
    my $rc = snap("--date", $D, @args);
    my $ok = isnt($rc, 0, $label);
    $ok = like($lastOut, $re, "message") && $ok;
    $ok = is(ls_data(), $before, "nothing created, latest unchanged") && $ok;
    diag($lastOut) unless $ok;
}

# starts snapshot.pl in the background; output goes to control/bg.out
sub snap_bg {
    my @c = snap_cmd(@_);
    my $pid = fork();
    die "fork: $!" unless defined $pid;
    if (! $pid) {
        open(STDIN, '<', '/dev/null');
        open(STDOUT, '>', "$ctl/bg.out");
        open(STDERR, '>&', \*STDOUT);
        exec(@c[2 .. $#c]) or POSIX::_exit(127);
    }
    return $pid;
}

sub wait_for {
    my ($code) = @_;
    for (1 .. 600) {
        return 1 if $code->();
        sleep(0.05);
    }
    return 0;
}

sub run_prune {
    my $cmd = join(" ", map { shq($_) } ($^X, $prune, "--root", $root, @_));
    $lastOut = `$cmd </dev/null 2>&1`;
    return $? >> 8;
}

sub latest { return readlink("$root/data/latest") // "(none)"; }

sub clear_log { unlink("$ctl/access.log"); }

sub requests {
    my ($method, $re) = @_;
    open(my $fh, '<', "$ctl/access.log") or return 0;
    my $n = grep { /^$method / && /$re/ } <$fh>;
    return $n;
}
sub struct_gets { return requests("GET", "data/structures|compatible"); }
sub struct_heads { return requests("HEAD", "data/structures|compatible"); }

# everything in data/ except the lock file
sub ls_data {
    opendir(my $dh, "$root/data") or return "(no data directory)";
    my $l = join(" ", sort grep { !/^\.\.?$/ && $_ ne '.snapshot.lock' } readdir($dh));
    return $l." latest->".latest();
}

sub partials {
    opendir(my $dh, "$root/data") or return ();
    return map { "$root/data/$_" } sort grep { /\.partial$/ } readdir($dh);
}

# every file, directory and link below a directory, with what identifies it
sub tree {
    my ($top) = @_;
    my %t;
    find({ no_chdir => 1, wanted => sub {
        my @st = lstat($_);
        (my $rel = $_) =~ s{^\Q$top\E/?}{};
        return if ($rel eq '' || $rel =~ /\.snapshot\.lock$/);
        $t{$rel} = (-l _) ? "link ".readlink($_) : (-d _) ? "dir" : "file $st[7] $st[9]";
    } }, $top);
    return \%t;
}

# snapshot that holds the file a member link points to
sub source_of {
    my ($snap, $path) = @_;
    my $l = readlink("$root/data/$snap/snapshot/$path") // return "(no link)";
    return "(dangling)" unless (-f "$root/data/$snap/snapshot/$path");
    return ($l =~ m{^\.\./\.\./\.\./\.\./(\d{6})/files/\Q$path\E$}) ? $1 : "(odd link $l)";
}
sub mtime_of { return (stat("$root/data/$_[0]/snapshot/$_[1]"))[9]; }

sub manifest {
    my ($snap) = @_;
    my %rows;
    open(my $fh, '<', "$root/data/$snap/files/MANIFEST.tsv") or return {};
    my @cols;
    while (<$fh>) {
        chomp;
        if (/^#path\t/) { s/^#//; @cols = split(/\t/); next; }
        next if /^#/;
        my %r;
        @r{@cols} = split(/\t/, $_, -1);
        $rows{$r{path}} = \%r;
    }
    return \%rows;
}

sub slurp {
    my ($f) = @_;
    open(my $fh, '<', $f) or return undef;
    binmode($fh);
    local $/;
    return scalar(<$fh>);
}

sub spew {
    my ($f, $bytes) = @_;
    open(my $fh, '>', $f) or die "$f: $!";
    binmode($fh);
    print $fh $bytes;
    close($fh);
}

sub gunzip_text {
    my ($f) = @_;
    my $q = shq($f);
    return scalar(`gzip -dc $q 2>/dev/null`);
}

# valid gzip data of exactly $len bytes whose content is not $text: a
# shorter text, padded with a file name in the gzip header
sub alt_gz {
    my ($text, $len) = @_;
    my @lines = split(/^/, $text);
    foreach my $drop (1 .. $#lines) {
        my $t = join("", @lines[0 .. $#lines - $drop]);
        my $pad = $len - length(gz($t)) - 1;
        next if ($pad < 1);
        my $out;
        gzip(\$t => \$out, Time => 0, Name => "x" x $pad) or die $GzipError;
        return $out if (length($out) == $len);
    }
    die "alt_gz: no data of $len bytes";
}

# replaces the content of a file, keeping its times
sub damage {
    my ($f, $bytes) = @_;
    my @st = stat($f) or die "$f: $!";
    spew($f, $bytes);
    utime($st[8], $st[9], $f) or die;
}

# checks a published snapshot as a whole
sub snapshot_ok {
    my ($snap, $label, %o) = @_;
    my $nStruct = $o{n} // nstruct();
    my $d = "$root/data/$snap";
    my $m = manifest($snap);
    my @struct = grep { $m->{$_}{kind} !~ /^(seq|index)$/ } keys %$m;
    my @bad;
    push(@bad, "no COMPLETE marker") unless ($o{old_style} || -f "$d/COMPLETE");
    push(@bad, scalar(@struct)." structure rows, expected $nStruct") if (@struct != $nStruct);
    foreach my $p (@struct) {
        my $r = $m->{$p};
        my $f = "$root/data/$r->{storage}";
        my $via = "$d/snapshot/$p";
        push(@bad, "$p: link does not resolve") unless (-f $via);
        (my $all = $p) =~ s{^([a-z-]*)hash/../}{${1}all/};
        push(@bad, "$p: 'all' link does not resolve") unless (-f "$d/snapshot/$all");
        push(@bad, "$p: link is not the manifest storage")
            unless ((readlink($via) // '') eq "../../../../$r->{storage}");
        my $bytes = slurp($f) // '';
        push(@bad, "$p: size") unless ($r->{size} == length($bytes));
        push(@bad, "$p: md5") unless ($r->{md5} eq '-' || $r->{md5} eq md5_hex($bytes));
        push(@bad, "$p: sha256")
            unless ($r->{sha256} eq '-' || $r->{sha256} eq sha256_hex(gunzip_text($f)));
        push(@bad, "$p: status/source") unless
            ($r->{source} eq substr($r->{storage}, 0, 6) &&
             $r->{status} eq ($r->{source} eq $snap ? 'fetched' : 'reused'));
    }
    # ParsePDB reads everything in files/hash and files/obs-hash
    foreach my $h ('hash', 'obs-hash') {
        opendir(my $dh, "$d/files/$h") or next;
        foreach my $n (grep { !/^\.\.?$/ } readdir($dh)) {
            push(@bad, "unexpected $h/$n") unless ($n eq 'removed' || $n =~ /^[a-z0-9]{2}$/);
        }
    }
    my @junk;
    find(sub { push(@junk, $File::Find::name) if /\.tmp\.|\.partial$|^HEADS\.tsv$|^latest\.tmp/ }, "$root/data");
    push(@bad, "temporary files left: @junk") if @junk;
    foreach my $n (@indexNames) {
        push(@bad, "index $n.gz") unless ((gunzip_text("$d/files/index/$n.gz") // '') eq "contents of $n\n");
    }
    foreach my $n ("pdbseqres", "pdbprot", "onhold") {
        push(@bad, "seq $n") unless (-s "$d/snapshot/seq/$n.$snap.fa" && -s "$d/snapshot/seq/$n.fa");
    }
    push(@bad, "holdings copy") unless ($o{old_style} || -s "$d/files/holdings/all_removed_entries.json.gz");
    ok(! @bad, $label) or diag(join("\n", @bad));
}

# ---------------------------------------------------------------- server

set_faults();
my $serverPid = fork();
die "fork: $!" unless defined $serverPid;
if (! $serverPid) {
    exec("python3", "$here/mockserver.py", $doc, $ctl) or POSIX::_exit(127);
}
END { kill('TERM', $serverPid) if $serverPid; }
for (1 .. 100) { last if (-s "$ctl/port"); sleep(0.1); }
my $port = slurp("$ctl/port") or die "mock server did not start\n";
$base = "http://127.0.0.1:$port";

my $old = '2029-12-14T10:00:00+00:00';
add_entry("1aa$_", date => $old) foreach ('a' .. 'k');
add_entry("1aal", date => $old, keys => [qw(mmcif pdbml pdb_bundle)]);
add_entry("1zza", date => $old, obsolete => 1);
add_entry("1zzb", date => $old, obsolete => 1);
# a removed entry with no date and its file somewhere else: not archived
add_entry("1zzc", obsolete => 1, keys => ['pdb'], date => undef);
$site{"1zzc"}{path}{pdb} = "data/structures/models/current/pdb/zz/pdb1zzc.ent.gz";
write_site();

# ---------------------------------------------------------------- cases

subtest "latest missing, no terminal, no --init: dies, no loop" => sub {
    next_date();
    my $rc = snap("--date", $D);
    isnt($rc, 0, "exit status $rc");
    isnt($rc, 124, "did not run into the timeout");
    like($lastOut, qr/--init/, "message names --init");
    ok(! -e "$root/data", "nothing created");
    ok(requests("GET", ".") == 0, "no request made");
};

subtest "--init with a failing guard creates nothing, not even the data directory" => sub {
    $hold{rename_pdb_key} = "pdb_legacy";
    write_site();
    isnt(snap("--init", "--date", $D), 0, "run fails");
    like($lastOut, qr/guard FAIL PDB-format share/, "guard reported");
    ok(! -e "$root/data", "no data directory, no lock file");
    %hold = ();
    write_site();
};

subtest "first run with --init (default guards)" => sub {
    publish("exit status", "--init");
    snapshot_ok($D, "snapshot complete and consistent");
    is(nstruct(), 42, "42 structure files");
    my $m = manifest($D);
    is(scalar(keys %$m), 42 + 3 + 11, "manifest rows (42 structure, 3 sequence, 11 index)");
    is(scalar(grep { $_->{status} ne 'fetched' || $_->{sha256} eq '-' } values %$m), 0,
       "all fetched, all with checksums");
    ok(! exists $m->{"obs-hash/zz/pdb1zzc.ent.gz"}, "file at an unexpected path not archived");
    like($lastOut, qr/pdb-obs\s+0\s+2\s+2\s+0\s+0\s+1\n/, "and counted as skipped");
    is(readlink("$root/hash"), "data/latest/snapshot/hash", "root links made");
    ok(-f "$root/all/pdb1aaa.ent.gz", "file reachable through root link");
    is(slurp("$root/seq/pdbprot.fa"), ">1aaa_A mol:protein length:6  PROTEIN ONE\nMKVLAA\n".
       ">1aab_A mol:protein length:5  PROTEIN TWO\nGGSAW\n", "protein sequences");
    is(slurp("$root/seq/onhold.fa"), $onhold, "on-hold sequences kept although no Last-Modified");
    ok(! partials(), "no staging directory left");
    ok(! glob("$tmp/*"), "temporary directory empty");
    unlike($lastOut, qr/pdb1aaa/, "no per-file output without --verbose");
};

subtest "rerun on the same day after success" => sub {
    refused("refused", qr/already exists/);
};

subtest "--init on an archive that has snapshots" => sub {
    unlink("$root/data/latest") or die;
    my $before = tree("$root/data");
    clear_log();
    isnt(snap("--init", "--date", snap_date($week + 1)), 0, "refused");
    like($lastOut, qr/--init given, but .* already holds 1 snapshot/, "message");
    is(requests("(GET|HEAD)", "."), 0, "no request made");
    is_deeply(tree("$root/data"), $before, "nothing changed");
    isnt(snap("--date", snap_date($week + 1)), 0, "and without --init it stops too");
    symlink($L, "$root/data/latest") or die;
};

subtest "second run, nothing changed: all reused (default guards)" => sub {
    next_date();
    clear_log();
    publish("exit status");
    is(struct_gets(), 0, "no structure file fetched");
    is(requests("HEAD", "."), 0, "no HEAD request");
    snapshot_ok($D, "snapshot consistent");
    my $m = manifest($D);
    is(scalar(grep { $_->{kind} !~ /seq|index/ && ($_->{status} ne 'reused' ||
                     $_->{source} ne snap_date(1) || $_->{sha256} eq '-') } values %$m), 0,
       "all rows reused from the first snapshot, checksums carried over");
    ok(! -e "$root/data/$D/files/hash/removed", "no removed list");
    like($lastOut, qr/total: fetched 0 .*reused 42, removed 0, failed 0/, "summary");
};

subtest "changed content, new entry, entry that vanishes from the holdings" => sub {
    $site{"1aab"}{files}{pdb} =~ s/version 1/version 2/;
    bump("1aab");
    $lm{rpath("1aab", "pdb")} = 1894000000;
    add_entry("1aam");
    delete $site{"1aac"};
    write_site();
    next_date();
    refused("default guards refuse an entry that is listed nowhere any more",
            qr/guard FAIL removed files whose entry is listed nowhere: 3 \(limit 0\)/);
    like($lastOut, qr/removed, entry listed nowhere: hash\/aa\/pdb1aac\.ent\.gz/, "and name the files");
    clear_log();
    publish("accepted with --max-vanished 3", "--max-vanished", 3);
    is(struct_gets(), 4, "fetched the changed file and the three new ones");
    is(struct_heads(), 3, "HEAD for the three files of the entry with a new date");
    is(source_of($D, "hash/aa/pdb1aab.ent.gz"), $D, "changed file is new");
    like(gunzip_text("$root/hash/aa/pdb1aab.ent.gz"), qr/^pdb file of 1aab, version 2\n/, "new content");
    is(mtime_of($D, "hash/aa/pdb1aab.ent.gz"), 1894000000, "modification time from Last-Modified");
    is(source_of($D, "cif-hash/aa/1aab.cif.gz"), snap_date(1), "its unchanged mmCIF is reused");
    is(source_of($D, "xml-hash/aa/1aam.xml.gz"), $D, "new entry fetched");
    is(slurp("$root/data/$D/files/hash/removed"), "pdb1aac.ent.gz\n", "removed list (PDB)");
    is(slurp("$root/data/$D/files/cif-hash/removed"), "1aac.cif.gz\n", "removed list (mmCIF)");
    like(slurp("$root/data/$D/files/REMOVED.tsv"),
         qr/^hash\/aa\/pdb1aac\.ent\.gz\t1aac\tpdb\t\d{6}\/files\/hash\/aa\/pdb1aac\.ent\.gz\t.*\tvanished$/m,
         "REMOVED.tsv records where the content is kept");
    ok(! -e "$root/data/$D/snapshot/hash/aa/pdb1aac.ent.gz", "removed entry not linked");
    snapshot_ok($D, "snapshot consistent");
};
my $removed1aac = $L;

subtest "Last-Modified changed, content identical: not fetched" => sub {
    bump("1aad");
    my $before = mtime_of($L, "hash/aa/pdb1aad.ent.gz");
    $lm{rpath("1aad", $_)} = 1895000000 foreach (qw(pdb mmcif pdbml));
    write_site();
    next_date();
    clear_log();
    publish("exit status");
    is(requests("HEAD", "1aad"), 3, "three HEAD requests");
    is(struct_gets(), 0, "nothing fetched");
    is(source_of($D, "hash/aa/pdb1aad.ent.gz"), snap_date(1), "still the file of the first snapshot");
    is(mtime_of($D, "hash/aa/pdb1aad.ent.gz"), $before, "modification time untouched");
    is(manifest($D)->{"hash/aa/pdb1aad.ent.gz"}{listed_date}, $site{"1aad"}{date},
       "new listed date recorded");
    snapshot_ok($D, "snapshot consistent");
};

subtest "multipart ETag: compared by content, then by the recorded ETag" => sub {
    set_faults(url_of("1aae", "mmcif") => { etag => "0123456789abcdef0123456789abcdef-2" });
    bump("1aae");
    write_site();
    next_date();
    clear_log();
    publish("exit status");
    is(requests("GET", "1aae.cif"), 1, "downloaded once to compare");
    is(source_of($D, "cif-hash/aa/1aae.cif.gz"), snap_date(1), "identical: old file kept");
    ok(! -e "$root/data/$D/files/cif-hash/aa/1aae.cif.gz", "download not kept");
    like($lastOut, qr/identical to the archived file, not kept: 1/, "reported");
    is(manifest($D)->{"cif-hash/aa/1aae.cif.gz"}{etag}, "0123456789abcdef0123456789abcdef-2",
       "ETag recorded");
    bump("1aae");
    write_site();
    next_date();
    clear_log();
    publish("next run");
    is(requests("HEAD", "1aae.cif"), 1, "HEAD");
    is(requests("GET", "1aae.cif"), 0, "same ETag as recorded: not downloaded");
    snapshot_ok($D, "snapshot consistent");
    set_faults();
};

subtest "multipart ETag that the local file reproduces: not downloaded" => sub {
    my $bytes = slurp("$pub/".rpath("1aaf", "mmcif"));
    set_faults(url_of("1aaf", "mmcif") => { etag => md5_hex(Digest::MD5::md5($bytes))."-1" });
    bump("1aaf");
    write_site();
    next_date();
    clear_log();
    publish("exit status");
    is(requests("HEAD", "1aaf"), 3, "three HEAD requests");
    is(struct_gets(), 0, "nothing downloaded");
    is(source_of($D, "cif-hash/aa/1aaf.cif.gz"), snap_date(1), "old file kept");
    like(manifest($D)->{"cif-hash/aa/1aaf.cif.gz"}{etag}, qr/^[0-9a-f]{32}-1$/, "ETag recorded");
    snapshot_ok($D, "snapshot consistent");
    set_faults();
};

subtest "same content compressed differently: not stored again" => sub {
    my $p = rpath("1aag", "pdb");
    my $before = slurp("$pub/$p");
    $gzLevel{$p} = 1;
    bump("1aag");
    write_site();
    isnt(slurp("$pub/$p"), $before, "the server's file has different bytes");
    my $mtime = mtime_of($L, "hash/aa/pdb1aag.ent.gz");
    next_date();
    clear_log();
    publish("exit status");
    is(requests("GET", "pdb1aag"), 1, "downloaded once to compare");
    is(source_of($D, "hash/aa/pdb1aag.ent.gz"), snap_date(1), "old file kept");
    is(mtime_of($D, "hash/aa/pdb1aag.ent.gz"), $mtime, "modification time untouched");
    ok(! -e "$root/data/$D/files/hash/aa/pdb1aag.ent.gz", "download not kept");
    like($lastOut, qr/same content compressed differently, not kept: 1/, "reported");
    like($lastOut, qr/total: fetched 0 /, "nothing counted as fetched");
    snapshot_ok($D, "snapshot consistent");
    delete $gzLevel{$p};
    write_site();
};

subtest "weak ETag, empty ETag, no Content-Length: content is compared" => sub {
    set_faults(url_of("1aah", "pdb") => { etag_raw => 'W/"'.md5_hex("x").'"' },
               url_of("1aai", "pdb") => { etag_raw => '""' },
               url_of("1aaj", "pdb") => { no_length => JSON::PP::true });
    bump("1aah", "1aai", "1aaj");
    write_site();
    next_date();
    refused("default guards refuse so many candidates", qr/guard FAIL pdb candidates 3 /);
    clear_log();
    publish("exit status", "--max-candidates", 0.5);
    is(requests("GET", "pdb1aa[hij]"), 3, "each of the three downloaded once");
    is(struct_gets(), 3, "and nothing else");
    is(scalar(grep { source_of($D, "hash/aa/pdb1aa$_.ent.gz") ne snap_date(1) } qw(h i j)), 0,
       "none stored again");
    like($lastOut, qr/identical to the archived file, not kept: 3/, "reported");
    unlike(manifest($D)->{"hash/aa/pdb1aah.ent.gz"}{etag}, qr{W/|"}, "a weak ETag is not recorded");
    snapshot_ok($D, "snapshot consistent");
    set_faults();
};

subtest "ETag of 32 hex digits that is not the MD5: --etag-not-md5" => sub {
    set_faults("*" => { etag => "00000000000000000000000000000000" });
    bump("1aak");
    write_site();
    next_date();
    clear_log();
    isnt(snap("--date", $D), 0, "run fails");
    like($lastOut, qr/differs from ETag 0{32} \(if the server's ETags are no longer MD5 sums, use --etag-not-md5\)/,
         "message names the option");
    ok(! -e "$root/data/$D", "nothing published");
    clear_log();
    publish("with --etag-not-md5", "--etag-not-md5");
    is(struct_gets(), 3, "the three candidates are downloaded");
    like($lastOut, qr/identical to the archived file, not kept: 3/, "and not stored: identical content");
    like($lastOut, qr/note: ETags were not used as MD5 sums/, "summary says so");
    snapshot_ok($D, "snapshot consistent");
    set_faults();
};

subtest "truncated body: detected, not kept; rerun next day" => sub {
    $site{"1aaf"}{files}{pdb} =~ s/version 1/version 2/;
    bump("1aaf");
    add_entry("1aax");
    write_site();
    set_faults(url_of("1aaf", "pdb") => { truncate => 20 });
    next_date();
    clear_log();
    my $before = latest();
    isnt(snap("--date", $D), 0, "run fails");
    like($lastOut, qr/FAILED hash\/aa\/pdb1aaf\.ent\.gz/, "failure reported");
    is(requests("GET", "pdb1aaf"), 3, "three attempts (two retries), then stop");
    ok(! -e "$root/data/$D", "no dated directory");
    is(latest(), $before, "latest unchanged");
    my ($stg) = partials();
    ok($stg && ! glob("$stg/files/hash/aa/*1aaf*"), "truncated file not kept");
    ok($stg && -f "$stg/files/hash/aa/pdb1aax.ent.gz", "other fetched files wait in the staging directory");
    # the next day the server is repaired, and 1aax is no longer listed
    set_faults();
    delete $site{"1aax"};
    write_site();
    next_date();
    publish("rerun succeeds");
    is(source_of($D, "hash/aa/pdb1aaf.ent.gz"), $D, "changed file fetched");
    ok(! glob("$root/data/$D/files/*/aa/*1aax*"), "files no longer listed are not in the snapshot");
    ok(! partials(), "no staging directory left");
    snapshot_ok($D, "snapshot consistent");
};

subtest "structure file without Last-Modified: the file fails" => sub {
    add_entry("1aan");
    write_site();
    set_faults(url_of("1aan", "pdb") => { no_last_modified => JSON::PP::true });
    next_date();
    my $before = latest();
    isnt(snap("--date", $D), 0, "run fails");
    like($lastOut, qr/FAILED hash\/aa\/pdb1aan\.ent\.gz: .*no usable Last-Modified/, "reported");
    ok(! -e "$root/data/$D", "no dated directory");
    is(latest(), $before, "latest unchanged");
    ok(! glob((partials())[0]."/files/hash/aa/*1aan*"), "file not kept");
    set_faults();
    publish("after the server sends it again");
    snapshot_ok($D, "snapshot consistent");
};

subtest "unparseable date: dies before creating anything" => sub {
    my $good = $site{"1aag"}{date};
    $site{"1aag"}{date} = "last Tuesday";
    write_site();
    next_date();
    refused("run fails", qr/Unparseable date 'last Tuesday'/);
    $site{"1aag"}{date} = $good;
    write_site();
};

subtest "guards (defaults): no pdb key, collapsed count, deletions" => sub {
    my $rm0 = scalar(grep { /removed|REMOVED/ } keys %{ tree("$root/data") });
    $hold{rename_pdb_key} = "pdb_legacy";
    write_site();
    refused("holdings without the pdb key: run fails", qr/guard FAIL PDB-format share 0\.0000/);
    like($lastOut, qr/guard FAIL removed files whose entry is still listed: 14 /, "still-listed guard");
    %hold = ();

    my %saved = %site;
    delete $site{"1aa$_"} foreach ('d' .. 'n');
    write_site();
    refused("collapsed entry count: run fails", qr/guard FAIL cif count 13 -> 2/);

    %site = %saved;
    delete $site{"1aa$_"} foreach ('i' .. 'k');
    write_site();
    refused("three entries removed, limit two: run fails", qr/guard FAIL pdb removed 3 \(limit 2\)/,
            "--max-deletions", 2, "--max-vanished", 99, "--max-count-change", 0.9);
    unlike($lastOut, qr/guard FAIL (?!pdb removed|cif removed|xml removed)/, "the only guard that fails");

    is(scalar(grep { /removed|REMOVED/ } keys %{ tree("$root/data") }), $rm0,
       "no removed list was written anywhere, staging included");
    %site = %saved;
    write_site();
};

subtest "entry loses its PDB-format file; another entry becomes obsolete" => sub {
    delete $site{"1aad"}{files}{pdb};
    $site{"1aae"}{obsolete} = 1;
    write_site();
    next_date();
    # the small test site falls below the default PDB-format share here
    my @o = ("--max-count-change", 0.6, "--min-pdb-share", 0.8);
    refused("default guards refuse", qr/guard FAIL removed files whose entry is still listed: 1 \(limit 0\)/, @o);
    like($lastOut, qr/removed, entry still listed \(not written to 'removed'\): hash\/aa\/pdb1aad\.ent\.gz/,
         "the file is named");
    publish("accepted with --max-still-listed 1", "--max-still-listed", 1, @o);
    is(slurp("$root/data/$D/files/hash/removed"), "pdb1aae.ent.gz\n",
       "'removed' has the obsoleted entry only, not the one that is still current");
    my $rm = slurp("$root/data/$D/files/REMOVED.tsv");
    like($rm, qr/^hash\/aa\/pdb1aad\.ent\.gz\t.*\tstill-listed$/m, "REMOVED.tsv has the lost file");
    like($rm, qr/^hash\/aa\/pdb1aae\.ent\.gz\t.*\tmoved$/m, "and the obsoleted one");
    like($lastOut, qr/note: removed but not written to 'hash\/removed' \(entry still listed\): pdb1aad\.ent\.gz/,
         "summary names it");
    ok(! -e "$root/data/$D/snapshot/hash/aa/pdb1aad.ent.gz", "lost file not in the snapshot");
    ok(-f "$root/data/$D/snapshot/obs-hash/aa/pdb1aae.ent.gz", "obsolete file archived");
    snapshot_ok($D, "snapshot consistent");
};

subtest "listed file returns 404: bounded failure, rerun the same day asks nothing twice" => sub {
    add_entry("1aao");
    add_entry("1aap");
    bump("1aaa");
    write_site();
    set_faults(url_of("1aap", "pdb") => { status => 404 });
    next_date();
    clear_log();
    my $before = latest();
    my @o = ("--max-count-change", 0.6, "--min-pdb-share", 0.8);
    isnt(snap("--date", $D, @o), 0, "run fails");
    is(requests("GET", "pdb1aap"), 3, "three attempts, then stop");
    is(requests("HEAD", "1aaa"), 3, "the candidate entry was checked");
    like($lastOut, qr/404/, "status reported");
    ok(! -e "$root/data/$D", "no dated directory");
    is(latest(), $before, "latest unchanged");
    is(scalar(partials()), 1, "staging directory kept for the rerun");
    set_faults();
    clear_log();
    publish("rerun succeeds without cleanup", @o);
    is(requests("GET", "1aao"), 0, "files fetched by the failed run are not fetched again");
    is(requests("HEAD", "."), 0, "and no HEAD request is repeated");
    is(struct_gets(), 1, "only the file that failed is fetched");
    like($lastOut, qr/kept from an interrupted run: 5/, "reported");
    ok(! partials(), "no staging directory left");
    snapshot_ok($D, "snapshot consistent");
};

subtest "--max-failures stops the run" => sub {
    add_entry("1ab$_") foreach ('a' .. 'f');
    write_site();
    set_faults(map { url_of("1ab$_", "pdb") => { status => 404 } } ('a' .. 'f'));
    next_date();
    clear_log();
    my @o = ("--max-count-change", 0.6);
    isnt(snap("--date", $D, @o, "--max-failures", 1), 0, "run fails");
    like($lastOut, qr/STOPPED: more than 1 files failed/, "stopped");
    is(scalar(() = $lastOut =~ /^FAILED /mg), 2, "after the second failure");
    is(requests("GET", "pdb1ab"), 6, "six requests: two files, three attempts each");
    ok(! -e "$root/data/$D", "nothing published");
    set_faults();
};

subtest "run killed half-way (SIGKILL), rerun the same day" => sub {
    my @o = ("--max-count-change", 0.6);
    set_faults("*" => { delay => 0.25 });
    my $before = latest();
    my $pid = snap_bg("--date", $D, @o, "--fresh");
    my $seen = wait_for(sub { my @f = glob("$root/data/.$D.partial/files/*/ab/*.gz"); @f >= 2 });
    kill('KILL', $pid);
    waitpid($pid, 0);
    ok($seen, "killed while fetching");
    ok(! -e "$root/data/$D", "no dated directory");
    is(latest(), $before, "latest unchanged");
    my @left = glob("$tmp/past.*");
    is(scalar(@left), 1, "the killed run could not remove its temporary directory");
    utime(time() - 3 * 86400, time() - 3 * 86400, @left);
    my @keptFiles = glob("$root/data/.$D.partial/files/*/ab/*.gz");
    my $kept = scalar(@keptFiles);
    set_faults();
    clear_log();
    publish("rerun succeeds without cleanup", @o);
    ok(! glob("$tmp/*"), "old temporary directory of the killed run removed");
    like($lastOut, qr/kept from an interrupted run: $kept\n/, "the $kept files of the killed run are kept");
    is(requests("GET", "1ab"), 18 - $kept, "exactly the other ".(18 - $kept)." are fetched");
    ok(! partials(), "no staging directory left");
    snapshot_ok($D, "snapshot consistent, no temporary files");
};

subtest "SIGTERM and SIGINT: clean stop, staging kept" => sub {
    add_entry("1aca");
    write_site();
    next_date();
    my $before = latest();
    foreach my $sig ('TERM', 'INT') {
        set_faults("*" => { delay => 0.3 });
        clear_log();
        my $pid = snap_bg("--date", $D);
        # in the middle of the sequence and index files
        my $seen = wait_for(sub { requests("GET", "derived_data/index/.") >= 2 });
        kill($sig, $pid);
        waitpid($pid, 0);
        my $rc = $? >> 8;
        ok($seen && $rc != 0, "SIG$sig: run ends with exit status $rc");
        like(slurp("$ctl/bg.out"), qr/Interrupted by signal $sig/, "SIG$sig: says why");
        my @junk;
        find(sub { push(@junk, $_) if /\.tmp/ }, "$root/data");
        ok(! @junk, "SIG$sig: no temporary file left in the archive") or diag("@junk");
        ok(! glob("$tmp/*"), "SIG$sig: temporary directory removed");
        is(scalar(partials()), 1, "SIG$sig: staging directory kept");
        is(latest(), $before, "SIG$sig: latest unchanged");
    }
    set_faults();
    publish("rerun succeeds");
    snapshot_ok($D, "snapshot consistent");
};

subtest "lock held by another run" => sub {
    open(my $fh, '>>', "$root/data/.snapshot.lock") or die;
    flock($fh, LOCK_EX) or die;
    next_date();
    refused("run refused", qr/holds the lock/);
    close($fh);
};

subtest "more than one staging directory" => sub {
    mkdir("$root/data/.300101.partial") or die;
    mkdir("$root/data/.300102.partial") or die;
    refused("run refused", qr/More than one interrupted run found/);
    rmdir("$root/data/.300101.partial") or die;
    rmdir("$root/data/.300102.partial") or die;
};

subtest "index listing lacks a required file; on-hold URL returns a web page" => sub {
    my $before = latest();
    bump("1aab");
    $indexHidden{"entries.idx"} = 1;
    write_site();
    clear_log();
    isnt(snap("--date", $D), 0, "run fails");
    like($lastOut, qr/lacks: entries\.idx/, "missing index file named");
    is(struct_heads() + struct_gets(), 0, "before any structure file is asked for");
    is(latest(), $before, "latest unchanged");
    %indexHidden = ();
    my $good = $onhold;
    $onhold = "<!DOCTYPE html><html><body>Moved</body></html>\n";
    write_site();
    isnt(snap("--date", $D), 0, "run fails");
    like($lastOut, qr/on-hold file is not FASTA/, "on-hold content check");
    is(latest(), $before, "latest unchanged");
    ok(! -e "$root/data/$D", "no dated directory");
    $onhold = $good;
    # a file that is already gzipped appears in the listing
    $indexExtra{"extra.idx.gz"} = gz("contents of extra.idx\n");
    write_site();
    publish("run with a new .gz name in the index listing");
    is(gunzip_text("$root/data/$D/files/index/extra.idx.gz"), "contents of extra.idx\n",
       "stored as it is, not compressed twice");
    snapshot_ok($D, "snapshot consistent");
    $indexExtra{"extra.idx"} = "plain\n";
    write_site();
    next_date();
    isnt(snap("--date", $D), 0, "listing with both a file and its .gz: run fails");
    like($lastOut, qr/both a file and its \.gz: extra\.idx\.gz/, "message");
    %indexExtra = ();
    write_site();
};

subtest "--dry-run creates nothing; --list names every file" => sub {
    add_entry("1acc");
    bump("1aah");
    delete $site{"1abf"};
    write_site();
    mkdir("$tmp/past.abcdefgh") or die;
    utime(time() - 3 * 86400, time() - 3 * 86400, "$tmp/past.abcdefgh");
    clear_log();
    my $before = tree($root);
    my $list = "$work/list.tsv";
    is(snap("--date", $D, "--dry-run", "--max-vanished", 3, "--list", $list), 0, "exit status") or diag($lastOut);
    like($lastOut, qr/^pdb\s+19\s+19\s+1\s+1\s+1\s+0$/m, "counts for PDB files: 1 new, 1 candidate, 1 removed");
    like($lastOut, qr/guard ok\s+pdb removed 1/, "guards reported");
    is_deeply(tree($root), $before, "archive unchanged: every name, size, modification time and link");
    is(requests("GET", "holdings"), 4, "four holdings files fetched");
    is(requests("(GET|HEAD)", "."), 4, "nothing else requested");
    ok(-d "$tmp/past.abcdefgh", "an old temporary directory is not touched without the lock");
    rmdir("$tmp/past.abcdefgh");
    ok(! glob("$tmp/*"), "its own temporary directory is removed");
    is(join("", sort grep { /pdb/ } split(/^/, slurp($list))),
       "candidate\thash/aa/pdb1aah.ent.gz\t\n".
       "new\thash/ac/pdb1acc.ent.gz\t\n".
       "removed\thash/ab/pdb1abf.ent.gz\tentry is listed nowhere; written to hash/removed\n",
       "--list names the new, candidate and removed PDB files");
    is(scalar(() = slurp($list) =~ /\n/g), 1 + 9, "nine names in all");
    is(snap("--date", $D, "--dry-run"), 2, "dry run with a failing guard exits 2");
    publish("real run", "--max-vanished", 3);
    snapshot_ok($D, "snapshot consistent");
};

subtest "entry that was removed is listed again: recognised, not stored again" => sub {
    add_entry("1aac", date => $old);
    delete $site{"1abe"};
    write_site();
    next_date();
    my @o = ("--max-vanished", 3);
    refused("default guards refuse", qr/guard FAIL new files that an earlier snapshot removed: 3 \(limit 0\)/, @o);
    like($lastOut, qr/reappeared \(removed in $removed1aac\): hash\/aa\/pdb1aac\.ent\.gz/, "the files are named");
    my $list = "$work/list2.tsv";
    snap("--date", $D, "--dry-run", "--list", $list, @o);
    like(slurp($list), qr/^reappeared\thash\/aa\/pdb1aac\.ent\.gz\tremoved in $removed1aac$/m, "also by --list");
    clear_log();
    publish("accepted with --max-reappeared 3", @o, "--max-reappeared", 3);
    is(requests("HEAD", "1aac"), 3, "three HEAD requests");
    is(struct_gets(), 0, "nothing fetched: the server has what the archive keeps");
    is(source_of($D, "hash/aa/pdb1aac.ent.gz"), snap_date(1), "linked to the copy of the first snapshot");
    ok(! glob("$root/data/$D/files/*/aa/*1aac*"), "nothing in files/ of this snapshot");
    like(slurp("$root/data/$D/files/REAPPEARED.tsv"),
         qr/^hash\/aa\/pdb1aac\.ent\.gz\t$removed1aac\treused \d{6}\/files\/hash\/aa\/pdb1aac\.ent\.gz$/m,
         "REAPPEARED.tsv lists it");
    like($lastOut, qr/hash\/aa\/pdb1aac\.ent\.gz: removed in $removed1aac, listed again; identical to .* ParsePDB will not see it/,
         "summary says what ParsePDB will not see");
    snapshot_ok($D, "snapshot consistent");
};
my $removed1abe = $L;

subtest "previous snapshot without a manifest; --manifest-only" => sub {
    # make the latest snapshot look like one of the old tool
    my $P = $L;
    unlink("$root/data/$P/files/MANIFEST.tsv", "$root/data/$P/files/REMOVED.tsv",
           "$root/data/$P/files/REAPPEARED.tsv", "$root/data/$P/COMPLETE");
    remove_tree("$root/data/$P/files/holdings");
    # entries listed as modified since 14 days before the previous snapshot
    my $cutoff = POSIX::strftime("%Y-%m-%d", gmtime($t0 + 7 * 86400 * $week - 14 * 86400));
    $site{"1aai"}{files}{pdb} =~ s/version 1/version 2/;
    bump("1aai");                                         # after the previous snapshot
    $site{"1aaj"}{date} = "${cutoff}T00:00:00+00:00";     # on the margin: checked
    my @inside = grep { defined $site{$_}{date} && $site{$_}{date} ge $cutoff } keys %site;
    my $nHead = 0;
    $nHead += scalar(keys %{$site{$_}{files}}) foreach (@inside);
    write_site();
    next_date();
    refused("default guards refuse so many candidates", qr/guard FAIL pdb candidates/);
    clear_log();
    publish("exit status", "--max-candidates", 0.9);
    like($lastOut, qr/membership from its link tree/, "membership read from the link tree");
    is(requests("HEAD", "1aai|1aaj"), 6, "HEAD for entries listed as modified since the margin");
    is(requests("HEAD", "1aak"), 0, "not for an entry listed as modified long before");
    is(requests("HEAD", "."), $nHead, "$nHead HEAD requests in all (".scalar(@inside)." entries inside the margin)");
    is(struct_gets(), 1, "only the changed file fetched");
    is(source_of($D, "hash/aa/pdb1aai.ent.gz"), $D, "changed file is new");
    snapshot_ok($D, "snapshot consistent");
    my $np = nstruct() - 1;
    my $m = manifest($D);
    is(scalar(grep { $_->{sha256} eq '-' } values %$m), $np, "reused files have no checksum yet");
    like($lastOut, qr/$np without checksums \(run: snapshot\.pl --manifest-only $D\)/, "and the summary says so");

    is(snap("--manifest-only", $D), 0, "--manifest-only completes it") or diag($lastOut);
    $m = manifest($D);
    is(scalar(grep { $_->{sha256} eq '-' || $_->{md5} eq '-' } values %$m), 0, "no row without checksums");
    snapshot_ok($D, "checksums are right");
    ok(! -e "$root/data/$D/files/MANIFEST.tsv.partial", "work file removed");

    is(snap("--manifest-only", $P, "--reference", "$root/data/$D/files/MANIFEST.tsv"), 0,
       "--manifest-only builds the manifest of the old-style snapshot") or diag($lastOut);
    like($lastOut, qr/ rows: 15 files read, $np checksums taken from other manifests/, "reusing known checksums");
    snapshot_ok($P, "and it is right", old_style => 1, n => $np + 1);
    unlink("$root/data/$P/files/MANIFEST.tsv");

    # a corrupt archived file is found
    my $Q = snap_date($week - 3);
    my $victim = "$root/data/".snap_date(1)."/files/cif-hash/aa/1aaa.cif.gz";
    my $bytes = slurp($victim);
    my @st = stat($victim);
    rename("$root/data/$Q/files/MANIFEST.tsv", "$work/saved-manifest") or die;
    spew($victim, substr($bytes, 0, length($bytes) - 6));
    isnt(snap("--manifest-only", $Q), 0, "--manifest-only fails on a truncated archived file");
    like($lastOut, qr/fails gzip test: \d{6}\/files\/cif-hash\/aa\/1aaa\.cif\.gz/, "and names it");
    ok(! -e "$root/data/$Q/files/MANIFEST.tsv", "no manifest written");
    spew($victim, $bytes);
    utime($st[8], $st[9], $victim);
    is(snap("--manifest-only", $Q), 0, "after repair it resumes and finishes") or diag($lastOut);
    like($lastOut, qr/ rows: 1 files read/, "reading only what was left");
    # path, storage, size, md5, sha256, status and source must be those the run wrote
    my $pick = sub { join("\n", sort map { join("\t", (split(/\t/, $_, -1))[0, 3, 4, 5, 6, 10, 11]) }
                          grep { !/^#/ } split(/\n/, $_[0])) };
    is($pick->(slurp("$root/data/$Q/files/MANIFEST.tsv")), $pick->(slurp("$work/saved-manifest")),
       "rebuilt manifest agrees with the one the run wrote");
};

subtest "entry removed by an old-style snapshot is listed again" => sub {
    add_entry("1abe", date => $old);
    write_site();
    next_date();
    clear_log();
    publish("exit status", "--max-reappeared", 3);
    like($lastOut, qr/reappeared \(removed in $removed1abe\): hash\/ab\/pdb1abe\.ent\.gz/, "recognised from its 'removed' list");
    is(struct_gets(), 0, "nothing fetched");
    isnt(source_of($D, "hash/ab/pdb1abe.ent.gz"), $D, "linked to the copy the archive keeps");
    snapshot_ok($D, "snapshot consistent");
};

subtest "archived file missing: fetched again within --max-repaired, never silently in bulk" => sub {
    my $first = snap_date(1);
    my $f = "$root/data/$first/files/hash/aa/pdb1aaa.ent.gz";
    my $mtime = (stat($f))[9];
    # the server's copy has a later Last-Modified by now
    $lm{rpath("1aaa", "pdb")} = 1897000000;
    write_site();
    rename($f, "$f.aside") or die;
    next_date();
    clear_log();
    publish("one missing file: run succeeds");
    is(struct_gets(), 1, "one file fetched");
    is(source_of($D, "hash/aa/pdb1aaa.ent.gz"), $D, "stored in this snapshot");
    is(mtime_of($D, "hash/aa/pdb1aaa.ent.gz"), $mtime,
       "with the modification time recorded for the lost file, not the server's new one");
    like($lastOut, qr/damaged or missing archived files fetched again: 1/, "summary counts it");
    like($lastOut, qr/hash\/aa\/pdb1aaa\.ent\.gz: archived file was missing, fetched again; content and modification time as recorded/,
         "and names it");
    snapshot_ok($D, "snapshot consistent");
    rename("$f.aside", $f) or die;

    # a whole directory of an old snapshot is gone
    my $dir = "$root/data/$first/files/cif-hash";
    rename($dir, "$dir.aside") or die;
    next_date();
    clear_log();
    my $before = latest();
    isnt(snap("--date", $D), 0, "many missing files: run fails");
    like($lastOut, qr/STOPPED: more than 5 archived files are missing/, "stopped at the limit");
    is(struct_gets(), 5, "after fetching five");
    ok(! -e "$root/data/$D", "nothing published");
    is(latest(), $before, "latest unchanged");
    rename("$dir.aside", $dir) or die;

    SKIP: {
        skip("running as root", 3) if ($> == 0);
        my $sub = "$root/data/$first/files/xml-hash/aa";
        chmod(0000, $sub) or die;
        clear_log();
        isnt(snap("--date", $D), 0, "unreadable directory: run fails");
        like($lastOut, qr/FAILED xml-hash\/aa\/\S+: cannot stat .*Permission denied/, "as an error, not as a missing file");
        is(struct_gets(), 0, "nothing fetched in its place");
        chmod(0755, $sub) or die;
    }
    publish("with the archive whole again");
    like($lastOut, qr/stale file\(s\) of an interrupted run removed/, "files fetched by the failed runs are dropped");
    is(scalar(grep { $_->{source} eq $D } grep { $_->{kind} !~ /seq|index/ } values %{ manifest($D) }), 0,
       "nothing stored again");
    snapshot_ok($D, "snapshot consistent");
};

subtest "archived file with the wrong size or empty: fails, --refetch compares and replaces" => sub {
    my $first = snap_date(1);
    # manifest of the latest snapshot says a wrong size
    my $mf = "$root/data/$L/files/MANIFEST.tsv";
    my $text = slurp($mf);
    my $pathA = "hash/aa/pdb1aab.ent.gz";
    ($text =~ s/^(\Q$pathA\E\t[^\t]+\t[^\t]+\t[^\t]+\t)(\d+)/$1.($2 + 1)/me) or die;
    spew($mf, $text);
    # and an archived file is empty
    my $pathB = "xml-hash/aa/1aab.xml.gz";
    my $fB = "$root/data/$first/files/$pathB";
    my @st = stat($fB);
    my $bytesB = slurp($fB);
    spew($fB, "");
    utime($st[8], $st[9], $fB);
    next_date();
    clear_log();
    my $before = latest();
    isnt(snap("--date", $D), 0, "run fails");
    like($lastOut, qr/FAILED \Q$pathA\E: the archived file \S+ is \d+ bytes, its manifest says \d+\..*--refetch \Q$pathA\E/,
         "wrong size: named, with the remedy");
    like($lastOut, qr/FAILED \Q$pathB\E: the archived file \S+ is empty\..*--refetch \Q$pathB\E/,
         "empty file: named, with the remedy");
    is(struct_gets(), 0, "nothing fetched");
    ok(! -e "$root/data/$D", "nothing published");
    is(latest(), $before, "latest unchanged");
    clear_log();
    publish("with --refetch for both", "--refetch", $pathA, "--refetch", $pathB);
    is(struct_gets(), 2, "both downloaded");
    isnt(source_of($D, $pathA), $D, "wrong size: the archived file agrees with the server and is kept");
    like($lastOut, qr/\Q$pathA\E: archived file agrees with the server, kept/, "and the summary says so");
    is(source_of($D, $pathB), $D, "empty file: replaced by a copy in this snapshot");
    is(mtime_of($D, $pathB), $st[9], "with the recorded modification time");
    is(slurp("$root/data/$D/snapshot/$pathB"), $bytesB, "and the right content");
    snapshot_ok($D, "snapshot consistent");
    spew($fB, $bytesB);
    utime($st[8], $st[9], $fB);
};

subtest "--max-changed stops a run that stores too many changed files" => sub {
    foreach my $id (qw(1aaa 1aab 1aaf)) {
        $site{$id}{files}{mmcif} =~ s/version \d/version 9/;
        bump($id);
    }
    write_site();
    next_date();
    my $before = latest();
    isnt(snap("--date", $D, "--max-candidates", 0.5), 0, "run fails");
    like($lastOut, qr/STOPPED: more than 2 cif files have changed content \(--max-changed 0\.1\)/, "stopped");
    ok(! -e "$root/data/$D", "nothing published");
    is(latest(), $before, "latest unchanged");
    clear_log();
    publish("accepted with --max-changed 0.5", "--max-candidates", 0.5, "--max-changed", 0.5);
    is(struct_gets(), 0, "the rerun fetches nothing again");
    is(requests("HEAD", "divided/pdb/|mmCIF"), 0, "repeats no HEAD");
    is(struct_heads(), 3, "and sends the three that the stopped run had not come to");
    snapshot_ok($D, "snapshot consistent");
};

subtest "--check-all sends a HEAD for every file" => sub {
    next_date();
    clear_log();
    publish("exit status", "--check-all");
    my $nh = nstruct();
    like($lastOut, qr/^--check-all: a HEAD request for each of the $nh files already archived\..*--max-compared/m,
         "the volume is stated before anything is asked for");
    is(struct_heads(), nstruct(), nstruct()." HEAD requests");
    is(struct_gets(), 0, "nothing fetched");
    snapshot_ok($D, "snapshot consistent");
};

subtest "corrupt gzip with a multipart ETag: detected; --fresh" => sub {
    add_entry("1acd");
    write_site();
    my $rp = rpath("1acd", "pdbml");
    put("$pub/$rp", "this is not gzip data\n");
    set_faults("/pub/pdb/$rp" => { etag => "0123456789abcdef0123456789abcdef-3" });
    next_date();
    my $before = latest();
    isnt(snap("--date", $D), 0, "run fails");
    like($lastOut, qr/FAILED xml-hash\/ac\/1acd\.xml\.gz: .*gzip test failed/, "gzip test reported");
    is(latest(), $before, "latest unchanged");
    my ($stg) = partials();
    ok($stg && ! glob("$stg/files/xml-hash/ac/*1acd*"), "corrupt file not kept");
    ok($stg && -f "$stg/files/hash/ac/pdb1acd.ent.gz", "staging directory holds the other files");
    set_faults();
    write_site();
    clear_log();
    publish("rerun with --fresh succeeds", "--fresh");
    like($lastOut, qr/removed .*\.$D\.partial \(--fresh\)/, "staging directory discarded");
    is(requests("GET", "1acd"), 3, "all three files fetched again");
    snapshot_ok($D, "snapshot consistent");
};

subtest "a finished snapshot newer than 'latest': run stops; --adopt" => sub {
    my $P = $L;
    next_date();
    publish("run");
    my $N = $D;
    # 'latest' is moved back, as after a failure between rename and switch
    unlink("$root/data/latest");
    symlink($P, "$root/data/latest") or die;
    foreach my $d ($N, snap_date($week + 1)) {
        clear_log();
        my $before = ls_data();
        isnt(snap("--date", $d), 0, "run for $d stops");
        like($lastOut, qr/Snapshot $N in .* is newer than 'latest' \($P\).*--adopt $N/s, "and names --adopt");
        is(ls_data(), $before, "nothing changed, 'latest' not moved");
        is(requests("(GET|HEAD)", "."), 0, "no request made");
    }
    rename("$root/data/$N/COMPLETE", "$work/COMPLETE") or die;
    isnt(snap("--adopt", $N), 0, "--adopt refuses a snapshot without the COMPLETE marker");
    rename("$work/COMPLETE", "$root/data/$N/COMPLETE") or die;
    rename("$root/data/$N/snapshot/cif-all", "$work/cif-all") or die;
    isnt(snap("--adopt", $N), 0, "--adopt refuses a snapshot without its link tree");
    like($lastOut, qr/no link tree \(cif-all\)/, "message");
    rename("$work/cif-all", "$root/data/$N/snapshot/cif-all") or die;
    isnt(snap("--adopt", $P), 0, "--adopt refuses a snapshot that is not newer");
    is(latest(), $P, "latest still unchanged");
    is(snap("--adopt", $N), 0, "--adopt of the complete snapshot") or diag($lastOut);
    is(latest(), $N, "latest switched");
    is(requests("(GET|HEAD)", "."), 0, "without any request");
};

subtest "a kind the previous snapshot had no files of" => sub {
    delete $site{"1aal"}{files}{pdb_bundle};
    write_site();
    next_date();
    publish("the only bundle goes away", "--max-still-listed", 1, "--max-count-change", 1);
    $site{"1aal"}{files}{pdb_bundle} = "a new bundle\n";
    write_site();
    next_date();
    refused("bundles appear again: default guards refuse", qr/guard FAIL bundle count 0 -> 1/,
            "--max-reappeared", 1);
    publish("accepted with --allow-new-kind", "--max-reappeared", 1, "--allow-new-kind");
    is(source_of($D, "bundle-hash/aa/1aal-pdb-bundle.tar.gz"), $D, "changed content: fetched");
    snapshot_ok($D, "snapshot consistent");
};

subtest "--hash-reused computes missing checksums in the run" => sub {
    unlink("$root/data/$L/files/MANIFEST.tsv");
    next_date();
    publish("exit status", "--hash-reused");
    my $m = manifest($D);
    is(scalar(grep { $_->{sha256} eq '-' || $_->{md5} eq '-' } values %$m), 0, "no row without checksums");
    snapshot_ok($D, "and they are right");
};

subtest "--max-hours stops a run that takes too long" => sub {
    add_entry("1ace");
    write_site();
    next_date();
    my $before = latest();
    isnt(snap("--date", $D, "--max-hours", 0.0000001), 0, "run fails");
    like($lastOut, qr/STOPPED: running longer than \S+ hours/, "stopped");
    ok(! -e "$root/data/$D", "nothing published");
    is(latest(), $before, "latest unchanged");
    publish("rerun without the limit");
    snapshot_ok($D, "snapshot consistent");
};

subtest "output to a closed pipe does not stop the run" => sub {
    next_date();
    my $cmd = join(" ", map { shq($_) } snap_cmd("--date", $D));
    system("$cmd </dev/null 2>&1 | head -c 0");
    is(latest(), $D, "snapshot published");
    $L = $D;
    snapshot_ok($D, "snapshot consistent");
};

subtest "root and temporary directory with a space in their names" => sub {
    my ($r2, $t2) = ("$work/root two", "$work/tmp two");
    make_path($t2);
    my @o = ("--root", $r2, "--tmpdir", $t2, "--min-pdb-share", 0.5);
    is(snap("--init", "--date", "300101", @o), 0, "first run") or diag($lastOut);
    is(readlink("$r2/data/latest"), "300101", "published");
    mkdir("$t2/past.abcdefgh") or die;
    utime(time() - 3 * 86400, time() - 3 * 86400, "$t2/past.abcdefgh");
    add_entry("1acf");
    write_site();
    set_faults(url_of("1acf", "pdb") => { status => 404 });
    isnt(snap("--date", "300108", @o), 0, "second run fails");
    ok(-d "$r2/data/.300108.partial", "staging directory kept");
    ok(! -e "$t2/past.abcdefgh", "old temporary directory removed");
    set_faults();
    is(snap("--date", "300109", @o), 0, "rerun the next day") or diag($lastOut);
    like($lastOut, qr/resuming: .*\.300108\.partial/, "found and resumed the staging directory");
    like($lastOut, qr/kept from an interrupted run: 2/, "its files are kept");
    is(readlink("$r2/data/latest"), "300109", "published");
    isnt(snap("--init", "--date", "300110", @o), 0, "--init refused here too");
    delete $site{"1acf"};
    write_site();
};


subtest "--refetch on a truncated archived file: the good download replaces it" => sub {
    my $path = "hash/aa/pdb1aah.ent.gz";
    my $row = manifest($L)->{$path};
    my $f = "$root/data/$row->{storage}";
    my $good = slurp($f);
    my @st = stat($f);
    is($row->{size}, length($good), "the manifest is right about the file");
    damage($f, substr($good, 0, length($good) - 6));
    next_date();
    clear_log();
    my $before = latest();
    isnt(snap("--date", $D), 0, "run fails");
    like($lastOut, qr/FAILED \Q$path\E: the archived file \S+ is \d+ bytes, its manifest says \d+\..*--refetch \Q$path\E/,
         "truncated file named, with the remedy");
    is(struct_gets(), 0, "nothing fetched");
    is(latest(), $before, "latest unchanged");
    clear_log();
    publish("with --refetch", "--refetch", $path);
    is(struct_gets(), 1, "downloaded once");
    is(source_of($D, $path), $D, "the download is stored in this snapshot");
    is(slurp("$root/data/$D/snapshot/$path"), $good, "the snapshot has the whole file");
    is(mtime_of($D, $path), $st[9], "with the recorded modification time");
    my $new = manifest($D)->{$path};
    is($new->{size}, length($good), "manifest: size of the whole file");
    is($new->{md5}, md5_hex($good), "manifest: its MD5");
    like($lastOut, qr/\Q$path\E: archived file had the wrong size, fetched again; content and modification time as recorded/,
         "summary says what was done");
    is(length(slurp($f)), length($good) - 6, "the damaged file itself is not touched");
    snapshot_ok($D, "snapshot consistent");
    damage($f, $good);
};

subtest "corrupt archived file of the recorded size: --verify finds it, --refetch replaces it" => sub {
    my $path = "cif-hash/aa/1aah.cif.gz";
    my $P = $L;
    my $row = manifest($P)->{$path};
    my $f = "$root/data/$row->{storage}";
    my $good = slurp($f);
    my @st = stat($f);
    # other content, same size, valid gzip, same modification time
    my $badBytes = alt_gz($site{"1aah"}{files}{mmcif}, length($good));
    damage($f, $badBytes);
    isnt(gunzip_text($f), gunzip_text("$pub/".rpath("1aah", "mmcif")), "archived file has other content");
    is(-s $f, $row->{size}, "and the recorded size");

    # found by --manifest-only --verify
    my $mfTime = (stat("$root/data/$P/files/MANIFEST.tsv"))[9];
    isnt(snap("--manifest-only", $P, "--verify"), 0, "--manifest-only --verify fails");
    like($lastOut, qr/checksum differs from manifest \(md5, sha256\): \Q$row->{storage}\E/, "and names the file");
    like($lastOut, qr/verify: \d+ files read and compared with their recorded checksums, 1 differ/, "one file differs");
    is((stat("$root/data/$P/files/MANIFEST.tsv"))[9], $mfTime, "manifest not touched");
    # run again: only the file that failed is read again
    isnt(snap("--manifest-only", $P, "--verify"), 0, "it fails again");
    my $rest = scalar(keys %{ manifest($P) }) - 1;
    like($lastOut, qr/verify: 1 files read and compared with their recorded checksums, 1 differ; $rest taken as read by an interrupted pass/,
         "reading only the file that differed; the other $rest were read by the first pass");
    unlink("$root/data/$P/files/MANIFEST.tsv.partial") or die;

    # an ordinary run does not read the file and reuses it; a run with --verify stops
    next_date();
    clear_log();
    my $before = latest();
    isnt(snap("--date", $D, "--verify"), 0, "run with --verify fails");
    like($lastOut, qr/^--verify: \d+ archived files \(\d+ bytes\) will be read/m, "the volume is stated first");
    like($lastOut, qr/FAILED \Q$path\E: the archived file \S+ does not have the checksums recorded in its manifest\..*--refetch \Q$path\E/,
         "the file is named, with the remedy");
    is(scalar(() = $lastOut =~ /^FAILED /mg), 1, "no other file fails");
    is(struct_gets(), 0, "nothing fetched");
    ok(! -e "$root/data/$D", "nothing published");
    is(latest(), $before, "latest unchanged");

    # found in passing when the file has to be compared with a download
    set_faults(url_of("1aah", "mmcif") => { etag_raw => 'W/"x"' });
    bump("1aah");
    write_site();
    isnt(snap("--date", $D), 0, "run that downloads the file to compare fails");
    like($lastOut, qr/FAILED \Q$path\E: the archived file \S+ does not have the checksum recorded in its manifest, and is not what the server has\..*--refetch \Q$path\E/,
         "named, with the remedy");
    set_faults();

    # --refetch replaces it, and reports a good file as good
    my $other = "xml-hash/aa/1aab.xml.gz";
    clear_log();
    publish("with --refetch", "--refetch", $path, "--refetch", $other, "--refetch", "hash/zz/pdb9zzz.ent.gz");
    is(requests("GET", "1aah.cif"), 1, "downloaded once");
    is(source_of($D, $path), $D, "the download is stored in this snapshot");
    is(slurp("$root/data/$D/snapshot/$path"), $good, "with the right content");
    is(mtime_of($D, $path), $st[9], "and the recorded modification time");
    like($lastOut, qr/\Q$path\E: archived file did not have its recorded checksum, fetched again; content and modification time as recorded/,
         "summary says what was done");
    is(slurp($f), $badBytes, "the damaged file itself is not touched");
    is(requests("GET", "1aab.xml"), 1, "--refetch of a good file downloads it");
    isnt(source_of($D, $other), $D, "and keeps the archived file");
    like($lastOut, qr/\Q$other\E: archived file agrees with the server, kept$/m, "summary says so");
    like($lastOut, qr/--refetch hash\/zz\/pdb9zzz\.ent\.gz did nothing/, "--refetch of a name that is not there is reported");
    snapshot_ok($D, "snapshot consistent");

    # the new snapshot verifies; its manifest is not rewritten
    $mfTime = (stat("$root/data/$D/files/MANIFEST.tsv"))[9];
    utime($mfTime - 100, $mfTime - 100, "$root/data/$D/files/MANIFEST.tsv");
    is(snap("--manifest-only", $D, "--verify"), 0, "--manifest-only --verify of the new snapshot") or diag($lastOut);
    my $nrows = scalar(keys %{ manifest($D) });
    like($lastOut, qr/verify: $nrows files read and compared with their recorded checksums, 0 differ/,
         "all $nrows files read and compared");
    is((stat("$root/data/$D/files/MANIFEST.tsv"))[9], $mfTime - 100, "manifest not rewritten");
    ok(! -e "$root/data/$D/files/MANIFEST.tsv.partial", "work file removed");
    isnt(snap("--manifest-only", $D, "--verify", "--reference", "$root/data/$P/files/MANIFEST.tsv"), 0,
         "--verify refuses --reference");
    next_date();
    publish("a run with --verify on the repaired archive", "--verify");
    my $n = nstruct();
    like($lastOut, qr/verified against the previous manifest \(--verify\): $n archived files; 0 have no recorded checksum/,
         "all $n reused files verified");
    snapshot_ok($D, "snapshot consistent");
    damage($f, $good);
};

subtest "first run that failed: --init again resumes, --init --fresh starts again; --adopt without 'latest'" => sub {
    my ($r3, $r4) = ("$work/root3", "$work/root4");
    my @o = ("--min-pdb-share", 0.5);
    my $n = nstruct();
    set_faults(url_of("1aaa", "pdb") => { status => 404 });
    foreach my $r ($r3, $r4) {
        isnt(snap("--init", "--date", "300101", @o, "--root", $r), 0, "first run fails");
        ok(-d "$r/data/.300101.partial" && ! -e "$r/data/latest", "staging directory left, no 'latest'");
    }
    set_faults();
    isnt(snap("--date", "300102", @o, "--root", $r3), 0, "without --init the rerun stops");
    like($lastOut, qr/--init \(also after a first run that failed\)/, "and names --init");
    clear_log();
    is(snap("--init", "--date", "300102", @o, "--root", $r3), 0, "--init again publishes") or diag($lastOut);
    like($lastOut, qr/resuming: .*\.300101\.partial/, "resuming the staging directory");
    like($lastOut, qr/kept from an interrupted run: @{[ $n - 1 ]}\n/, "the ".($n - 1)." files of the failed run are kept");
    is(struct_gets(), 1, "only the file that failed is fetched");
    is(readlink("$r3/data/latest"), "300102", "latest made");
    is(readlink("$r3/hash"), "data/latest/snapshot/hash", "root links made");
    ok(! -e "$r3/data/.300101.partial" && ! -e "$r3/data/.300102.partial", "no staging directory left");

    clear_log();
    is(snap("--init", "--fresh", "--date", "300102", @o, "--root", $r4), 0, "--init --fresh publishes") or diag($lastOut);
    like($lastOut, qr/removed .*\.300101\.partial \(--fresh\)/, "staging directory discarded");
    is(struct_gets(), $n, "all $n files fetched");
    is(readlink("$r4/data/latest"), "300102", "latest made");

    # the first run published, but 'latest' and the root links were not made
    unlink("$r3/data/latest") or die;
    unlink("$r3/$_") foreach (qw(hash all seq index));
    my $state = sub { opendir(my $dh, "$r3/data") or die; join(" ", sort grep { !/^\./ } readdir($dh)) };
    my $before = $state->();
    isnt(snap("--init", "--date", "300103", @o, "--root", $r3), 0, "--init is refused: there is a snapshot");
    like($lastOut, qr/already holds 1 snapshot.*snapshot\.pl --adopt 300102/s, "and names --adopt");
    is($state->(), $before, "nothing created");
    isnt(snap("--adopt", "300101", "--root", $r3), 0, "--adopt of a name that is not a snapshot fails");
    mkdir("$r3/data/300105") or die;
    isnt(snap("--adopt", "300102", "--root", $r3), 0, "--adopt refuses when a newer dated directory exists");
    like($lastOut, qr/not the newest snapshot/, "message");
    rmdir("$r3/data/300105") or die;
    ok(! -l "$r3/data/latest", "still no 'latest'");
    is(snap("--adopt", "300102", "--root", $r3), 0, "--adopt of the published snapshot") or diag($lastOut);
    is(readlink("$r3/data/latest"), "300102", "latest made");
    ok(-f "$r3/hash/aa/pdb1aaa.ent.gz" && -s "$r3/seq/pdbprot.fa", "root links made");
    is(snap("--date", "300109", @o, "--root", $r3), 0, "and the next run works") or diag($lastOut);
};

subtest "staged file whose ETag is not an MD5: kept on a rerun only if it can be verified" => sub {
    add_entry("1adb");
    add_entry("1adc");
    write_site();
    my %f = (url_of("1adb", "mmcif") => { etag => "0123456789abcdef0123456789abcdef-2" },
             url_of("1adb", "pdbml") => { etag => "fedcba9876543210fedcba9876543210-2" },
             url_of("1adb", "pdb")   => { etag_raw => 'W/"abc"' });
    set_faults(%f, url_of("1adc", "pdb") => { status => 404 });
    next_date();
    isnt(snap("--date", $D), 0, "run fails");
    my ($stg) = partials();
    my $cif = "$stg/files/cif-hash/ad/1adb.cif.gz";
    ok(-f $cif && -f "$stg/files/xml-hash/ad/1adb.xml.gz" && -f "$stg/files/hash/ad/pdb1adb.ent.gz",
       "the three files of 1adb wait in the staging directory");
    # one of them is no longer what was fetched: same length, valid gzip
    my $server = slurp("$pub/".rpath("1adb", "mmcif"));
    damage($cif, alt_gz($site{"1adb"}{files}{mmcif}, length($server)));
    set_faults(%f);
    clear_log();
    publish("rerun");
    is(requests("HEAD", "."), 0, "no HEAD request is repeated");
    is(requests("GET", "1adb.cif"), 1, "multipart ETag, file changed since it was fetched: fetched again");
    is(slurp("$root/data/$D/snapshot/cif-hash/ad/1adb.cif.gz"), $server, "the snapshot has the server's file");
    is(requests("GET", "1adb.xml"), 0, "multipart ETag as when fetched, file unchanged: kept");
    is(requests("GET", "pdb1adb"), 1, "weak ETag: fetched again");
    is(requests("GET", "1adc"), 1, "MD5 ETags: kept; only the file that failed is fetched");
    like($lastOut, qr/kept from an interrupted run: 3\n/, "three files kept");
    snapshot_ok($D, "snapshot consistent");
    set_faults();
};

subtest "download with a multipart ETag: checked when it has one part" => sub {
    add_entry("1add");
    write_site();
    my $bytes = slurp("$pub/".rpath("1add", "pdbml"));
    set_faults(url_of("1add", "pdbml") => { etag => md5_hex(Digest::MD5::md5("other"))."-1" },
               url_of("1add", "mmcif") => { etag => "0123456789abcdef0123456789abcdef-2" });
    next_date();
    my $before = latest();
    isnt(snap("--date", $D), 0, "run fails");
    like($lastOut, qr/FAILED xml-hash\/ad\/1add\.xml\.gz: .*does not give the one-part ETag \S+-1 .*--etag-not-md5/,
         "a one-part ETag that the file does not give fails the download");
    like($lastOut, qr/multipart ETag that could not be checked against it .*: 1\n/,
         "a two-part ETag cannot be checked: accepted, and the summary counts it");
    is(latest(), $before, "latest unchanged");
    set_faults(url_of("1add", "pdbml") => { etag => md5_hex(Digest::MD5::md5($bytes))."-1" },
               url_of("1add", "mmcif") => { etag => "0123456789abcdef0123456789abcdef-2" });
    clear_log();
    publish("with the right one-part ETag");
    is(struct_gets(), 1, "only that file is fetched");
    is(source_of($D, "xml-hash/ad/1add.xml.gz"), $D, "and stored");
    snapshot_ok($D, "snapshot consistent");
    set_faults();
};

subtest "REMOVED.tsv with a bad storage column: the run stops" => sub {
    my ($with) = grep { -f "$root/data/$_/files/REMOVED.tsv" } reverse map { snap_date($_) } ($week - 20 .. $week);
    ok($with, "a recent snapshot has a REMOVED.tsv") or return;
    my $rf = "$root/data/$with/files/REMOVED.tsv";
    my $text = slurp($rf);
    (my $bad = $text) =~ s{^([^#\t][^\t]*\t[^\t]*\t[^\t]*\t)[^\t]+}{$1../../../etc/passwd}m or die;
    damage($rf, $bad);
    add_entry("1ade");
    write_site();
    next_date();
    refused("run refused", qr/REMOVED\.tsv line \d+: unexpected storage path/);
    damage($rf, $text);
    publish("with the file as it was");
    snapshot_ok($D, "snapshot consistent");
};

subtest "--max-compared stops a run that downloads files only to find them unchanged" => sub {
    bump("1aak");
    write_site();
    next_date();
    my $before = latest();
    clear_log();
    isnt(snap("--date", $D, "--etag-not-md5", "--max-compared", 1), 0, "run fails");
    like($lastOut, qr/STOPPED: more than 1 files were downloaded only to find that they are the archived copy \(--max-compared\)/,
         "stopped");
    is(struct_gets(), 2, "after the second such download");
    ok(! -e "$root/data/$D", "nothing published");
    is(latest(), $before, "latest unchanged");
    clear_log();
    publish("with the default limit", "--etag-not-md5");
    is(struct_gets(), 3, "the downloads that were not kept are made again");
    is(struct_heads(), 1, "the HEAD requests are not");
    snapshot_ok($D, "snapshot consistent");
};

subtest "makelinks.pl remakes a pruned link tree exactly" => sub {
    my $P = snap_date($week - 1);
    my $before = tree("$root/data/$P/snapshot");
    cmp_ok(scalar(keys %$before), '>', 150, "the tree of $P has ".scalar(keys %$before)." links and directories");
    is(run_prune($P), 0, "pruned") or diag($lastOut);
    ok(! -e "$root/data/$P/snapshot", "link tree gone");
    my $cmd = "PAST_ROOT=".shq($root)." ".shq($^X)." ".shq($makelinks)." $P";
    my $out = `$cmd 2>&1 </dev/null`;
    is($? >> 8, 0, "makelinks.pl $P") or diag($out);
    is_deeply(tree("$root/data/$P/snapshot"), $before, "identical tree: every link and its target");
    $out = `$cmd 2>&1 </dev/null`;
    isnt($? >> 8, 0, "refuses a snapshot that is not pruned");
};

subtest "frozen (read-only) snapshots are never changed" => sub {
    # freeze the latest snapshot the way release freezes are frozen
    my $F = $L;
    my $fdir = "$root/data/$F";
    system("chmod", "-R", "a-w", $fdir) == 0 or die;
    my $before = tree($fdir);
    my $thaw = sub { system("chmod", "-R", "u+w", $fdir); };

    $site{"1aai"}{files}{pdb} =~ s/version \d+/version 77/;
    bump("1aai");
    write_site();
    next_date();
    my $ok = publish("a run after a frozen snapshot publishes");
    is_deeply(tree($fdir), $before, "frozen snapshot unchanged (names, sizes, times, links)");
    is(source_of($D, "hash/aa/pdb1aac.ent.gz") eq $D ? "new" : "linked", "linked",
       "unchanged file is linked into an earlier snapshot");
    snapshot_ok($D, "new snapshot consistent") if $ok;

    isnt(run_prune($F), 0, "prune.pl refuses the frozen snapshot");
    like($lastOut, qr/frozen/, "and says why");
    isnt(snap("--manifest-only", $F), 0, "--manifest-only refuses to write into it");
    like($lastOut, qr/frozen/, "and says why");
    isnt(snap("--manifest-only", $F, "--verify"), 0, "so does --manifest-only --verify");
    isnt(snap("--manifest-only", $F, "--manifest-out", "$fdir/files/M2.tsv"), 0,
         "and --manifest-out inside the snapshot");
    is_deeply(tree($fdir), $before, "frozen snapshot still unchanged");

    # the refusal does not depend on the permission check failing: only the
    # top directory read-only, everything below writable
    $thaw->();
    chmod(0555, $fdir) or die;
    isnt(run_prune($F), 0, "prune.pl refuses when only the top directory is read-only");
    ok(-d "$fdir/snapshot/hash", "link tree still there");
    isnt(snap("--manifest-only", $F), 0, "--manifest-only refuses then too");
    like($lastOut, qr/frozen/, "and says why");
    my $out = "$tmp/frozen-manifest.tsv";
    unlink($out, "$out.partial");
    is(snap("--manifest-only", $F, "--verify", "--manifest-out", $out), 0,
       "--manifest-only --verify with --manifest-out elsewhere reads it") or diag($lastOut);
    my $cmd = "PAST_ROOT=".shq($root)." ".shq($^X)." ".shq($makelinks)." $F";
    my $mo = `$cmd 2>&1`;
    isnt($? >> 8, 0, "makelinks.pl refuses it");
    chmod(0755, $fdir) or die;
    unlink($out, "$out.partial");
};

subtest "prune.pl" => sub {
    my $l = latest();
    my $first = snap_date(1);
    my $before = ls_data();
    foreach my $arg ("latest", "$l", "$l/", "./$l", "..", ".", "$first/", "./$first",
                     "$root/data/$first", "$first/../$l", "", "30010", "3001011", "$first\n") {
        (my $show = $arg) =~ s/\n/\\n/;
        isnt(run_prune($arg), 0, "refuses '$show'");
    }
    symlink($l, "$root/data/999999") or die;
    isnt(run_prune("999999"), 0, "refuses a six-digit symlink to the latest snapshot");
    unlink("$root/data/999999");
    isnt(run_prune("300102"), 0, "refuses a snapshot that does not exist");
    ok(-d "$root/data/$l/snapshot/hash" && -d "$root/data/$first/snapshot/hash", "nothing was pruned");
    is(ls_data(), $before, "data directory unchanged");
    is(run_prune($first), 0, "prunes an old snapshot") or diag($lastOut);
    ok(! -e "$root/data/$first/snapshot", "link tree gone");
    ok(-f "$root/data/$first/files/hash/aa/pdb1aac.ent.gz", "files kept");
    ok(-f "$root/hash/aa/pdb1aac.ent.gz", "latest still reaches the file");
    isnt(run_prune($first), 0, "refuses a snapshot that is already pruned");
};

done_testing();
