PAST (PDB Archival Snapshot Toolkit) package 
--------------------------------------------

Copyright (C) 2001-2026 Degui Zhi and John-Marc Chandonia

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or (at
your option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
02110-1301, USA.

This file is part of the PAST (PDB Archival Snapshot Toolkit) package.

Contact:
scope@compbio.berkeley.edu

------------
Introduction
------------

Our PDB Archival Snapshot Toolkit (PAST) software allows users to
create date-stamped mirrors of the current PDB repository, and it
maintains all previous archives using an efficient storage mechanism.
Although the PDB development team now provides archival copies of the
PDB for selected dates, the archives are infrequent (updated yearly)
and do not provide an efficient means of storing the redundant data.
PAST is based on an efficient indexing technique similar to the ones
used by Picasa and iTunes.  We only store one physical copy for each
unique file.  PAST archival snapshots precisely mirror the directory
structure and files in the PDB, but all files in the snapshot
directory are replaced by links to a minimal set of data to prevent
wasteful use of space. PAST has been tested and deployed in our lab to
create weekly PDB snapshots corresponding to ASTEROIDS releases.

Important: If you are upgrading from PAST 1.21 or older, you need to
run "fixcompress.pl" on each of your archival snapshot directories,
starting with the earliest, in order to change from the PDB's old
compression standard (plain text or .Z) to the new standard (gzip).

-------------------------------
PAST releases
-------------------------------

PAST 2.0 (November 2026) makes a snapshot either complete and verified
  or not visible at all.  A file is fetched again only when its content
  changed (server MD5 and length, not dates); downloads are verified;
  every snapshot has a manifest with checksums; a failed run can simply
  be run again.  The BLAST databases are no longer built.  See "How
  snapshot.pl works" below.

PAST 1.42 was released in July 2021, and updates the sequence
  search URL at the PDB.

PAST 1.41 was released in March 2015, and contains a change to
  detect and work around truncated index files on the PDB's FTP site.

PAST 1.4 was released in January 2015, and contains changes to
  mirror the PDB's "PDB format-like files" (tar.gz format), which
  are released instead of PDB files for large structures.

PAST 1.31 was released in January 2013, and contains changes to
  be compatible with the wwPDB's new FTP server.

PAST 1.3 was released in October 2007, and contains changes to
  handle the wwPDB's new compression standards (e.g., all files gzipped),
  as well as better error recovery.

PAST 1.21 was released in July 2007, and fixed a bug in handling
  changes to the URL for the PDB's on-hold files.

PAST 1.2 was released in June 2007 to add the handling of PDBML files.

PAST 1.01 was released in June 2006 to handle changes to the PDB
  directory structure.

PAST 1.0 was released to the public in January 2006.

-------------------------------
Description of the PAST package
-------------------------------

This package includes programs for taking and maintaining archival
snapshots for the PDB. Currently it includes four perl programs:

1. snapshot.pl

This is the main program in this package. It makes a snapshot of the
PDB with today's date.

2. prune.pl

This program will prune back an old snapshot of the PDB to save on
inodes.  The "snapshot" directory is deleted, but the "files" will
still be there.  All info can be reconstructed with makelinks.pl.
The argument must be the six-digit name of a snapshot and nothing
else (no path, no trailing slash, not "latest"); the snapshot that
"latest" points to is never pruned.

3. makelinks.pl

This program will remake links from a pruned directory.  A snapshot
that has a manifest is remade from the manifest alone.  Otherwise the
previous snapshot must not have been pruned.  If the previous snapshot was
pruned, you have to run this program on it first (going all the way
back to the last un-pruned snapshot as necessary).

4. fixcompress.pl

This program is only needed when upgrading from PAST 1.21 or earlier
to PAST 1.3 or later.  You need to run it on each of your old archival
snapshots (starting with the earliest) in order to re-compress files
using gzip, which is now how the wwPDB stores its files.


------------
Requirements
------------

Linux, Perl 5.10 or later (the scripts start with #!/lab/bin/perl),
and gzip in the path.

snapshot.pl needs these Perl modules: LWP::UserAgent with
LWP::Protocol::https, HTTP::Date, Digest::MD5, Digest::SHA, File::Path,
File::Temp, and one JSON module (Cpanel::JSON::XS, JSON::XS,
JSON::MaybeXS or JSON::PP; JSON::PP is very slow on the holdings
file).  prune.pl and makelinks.pl need only modules that come with
Perl.  chromium and formatdb are no longer needed.

Memory: reading the holdings file takes about 1.5 GB (measured, with
a previous manifest of 790,000 rows loaded as well); allow 2 GB for a
whole run.


---------------
Getting started
---------------

First snapshot in a new archive:

    snapshot.pl --root /path/to/pdb --init

Every later snapshot:

    snapshot.pl

snapshot.pl never asks a question.  Without --init it stops if
<root>/data/latest is not a symlink to a snapshot; with --init it
stops if <root>/data/latest exists or <root>/data already holds a
snapshot (a dated directory).  A first run that failed leaves only
its staging directory: run again with --init, the same day or later,
and it is resumed like any other failed run (--init --fresh discards
it).  If a first run published its snapshot but could not make
"latest", "snapshot.pl --adopt <YYMMDD>" makes the link and the links
in <root>.

To see what a run would do, without fetching any structure file and
without creating anything in the archive:

    snapshot.pl --dry-run --list FILE

FILE gets one line per file that is new, a candidate for a check,
removed (with what will be done about it) or reappearing.


----------------------
How snapshot.pl works
----------------------

Two rules are kept throughout.  Content that is already in the
archive is never stored again and never gets a new modification time,
because ParsePDB makes a new release from either.  And anything that
would fetch again, or list as removed, more files than the limits
allow stops the run before a snapshot is published; every limit has
an option that raises it, to be given after looking at what the run
reported.

Layout.  <root>/data/<YYMMDD>/files holds the files first fetched for
that snapshot; <root>/data/<YYMMDD>/snapshot is a tree of links to
the files of this and earlier snapshots, which together are the PDB
as of that date.  <root>/data/latest points to the newest snapshot.

A run does this, in order:

1. Takes the lock <root>/data/.snapshot.lock, before it looks at
   anything in the archive.  A second run stops at once.  (The lock
   file stays; only the lock on it matters.)  --adopt takes the same
   lock.  Two exceptions: a dry run takes no lock, because it writes
   nothing; and the first run in a new archive, where <root>/data
   does not exist yet, makes the directory and takes the lock only
   after the guards have passed, and then looks again that nothing
   has appeared in it.

2. Stops if <root>/data holds a dated directory newer than "latest"
   (a snapshot that was finished but not switched to, or "latest" was
   moved back on purpose).  "snapshot.pl --adopt <YYMMDD>" makes such
   a snapshot the latest, if it is complete; otherwise move it away.

3. Fetches the four wwPDB holdings files and reads the membership of
   the previous snapshot: from its manifest if it has one, otherwise
   from one read of each "all" directory of its link tree.

4. Checks the guards.  If one fails, the run stops before anything is
   created, and "latest" is unchanged.  Per kind of file (pdb,
   pdb-obs, cif, cif-obs, xml, xml-obs, bundle):
     - a listed date that cannot be parsed, in the holdings or in the
       previous manifest;
     - the number of files changed by more than --max-count-change
       (0.10) of the previous number; a kind that had no files gets
       some (--allow-new-kind);
     - more than --max-deletions (60) files disappeared;
     - more than --max-candidates (0.20) of the listed files need a
       check (not with --check-all);
     - fewer than --min-pdb-share (0.90) of the current entries have
       a PDB-format file;
   and over all kinds:
     - a removed file whose entry is still in the same holdings file:
       the entry lost that format, or the file moved
       (--max-still-listed, 0);
     - a removed file whose entry is in neither holdings file
       (--max-vanished, 0);
     - a new file that one of the last --lookback (26) snapshots
       listed as removed (--max-reappeared, 0).
   The files of the last three guards are printed by name.

5. Works in <root>/data/.<YYMMDD>.partial.  First the four holdings
   files just fetched are copied to files/holdings.  Then the
   sequence file
   (compressed; stored uncompressed as seq/pdbseqres.<YYMMDD>.fa),
   the on-hold sequences (must be FASTA), seq/pdbprot.<YYMMDD>.fa,
   and every file in the listing of derived_data/index (stored
   gzipped; one that is listed gzipped is kept as it is).  The
   listing must contain all names in --index-required.

6. Then every listed structure file:
     - in the previous snapshot, and the date listed for its entry is
       the one recorded in the previous manifest: reused, no request;
     - listed date changed: a HEAD request.  The file is unchanged if
       the server's length equals the size of the archived file on
       disk and the server's ETag (the MD5 of the file) equals the
       MD5 recorded for it in the previous manifest, and then it is
       reused.  The archived file is not read for this; only when
       the manifest has no MD5 for it is the MD5 computed from the
       file.  Otherwise the file is downloaded and compared with the
       archived file as it is on disk (which is read then): it is
       stored only if neither the bytes nor the decompressed content
       are the same.  Last-Modified decides nothing.
     - not in the previous snapshot: fetched; but a file that an
       earlier snapshot listed as removed is first compared, in the
       same way, with the copy the archive still has, and linked to
       that copy if it is the same.
   An ETag of a multipart upload (32 hex digits, "-", a number; files
   over about 8 MB) is not an MD5.  Such a file is unchanged if its
   length and ETag equal those recorded in the previous manifest, or
   if the local file gives the same ETag when it is computed the way
   the server does (MD5 of the MD5s of 8 MiB parts); otherwise it is
   downloaded and compared.  A weak or empty ETag, or
   a missing Content-Length, also leads to download and comparison.
   --etag-not-md5 treats every ETag that way, for the day the
   server's ETags stop being MD5 sums (a download whose MD5 differs
   from its ETag fails otherwise).
   When the previous snapshot has no manifest (snapshots of PAST 1.x),
   a HEAD request is sent for every entry whose listed date is not
   earlier than --date-margin (14) days before the previous snapshot.
   --check-all sends a HEAD request for every file, and says how
   many before it starts.
   A reused file is therefore not read in an ordinary run: that it
   exists and has the recorded size is all that is checked.  --verify
   reads every reused file and compares its MD5 and the SHA-256 of
   its content with the previous manifest (it says first how many
   files and bytes that is); a file that differs fails, see "Damaged
   archive".  Files without recorded checksums are counted, not
   verified.

7. Every download goes to a temporary name and is renamed only after
   its length equals Content-Length, its MD5 equals the ETag (when
   the ETag is an MD5), gzip can decompress it, its modification time
   is set from Last-Modified, and it is synced to disk.  A structure
   file without Last-Modified fails, because ParsePDB reads the
   modification time.  A failed request is retried --retries (5)
   times, waiting --retry-delay (10) seconds and twice as long each
   time; the delay starts again for each file.  After a file has
   failed the run goes on with the others.  The run ends without
   publishing after --max-failures (20) failed files, when more than
   --max-changed (0.10 of the previous number, per kind) files were
   stored with changed content, when more than --max-repaired (5)
   archived files were missing, when more than --max-compared (100)
   files were downloaded only to turn out to be the archived copy
   (same bytes, or same content compressed differently), or after
   --max-hours (144).  The last limit is what stops a run from
   downloading the whole archive when the server's ETags stop being
   usable.
   A download with a multipart ETag is checked against it when it has
   one part ("-1"), where the part size does not matter, and fails if
   it differs.  With more parts it is checked for the usual part
   sizes; if none gives the ETag the download is accepted on its
   length and the gzip test, and the summary counts such downloads.

8. If nothing failed: removes from the staging directory the files
   that are not part of this snapshot; writes the "removed" lists,
   files/REMOVED.tsv, the manifest (files/MANIFEST.tsv) and the
   marker COMPLETE; renames the directory to <YYMMDD>; switches
   "latest" with one rename.
   What is synced to disk: every downloaded file before it gets its
   name; the copies of the holdings files, the sequence and index
   files, the lists, the manifest and the marker when they are
   written; and, before the rename, every directory of the staging
   directory that the run wrote into (the files/<kind>/<hh> and
   snapshot/<kind>/<hh> directories it used, the "all" directories
   and those above them), then <root>/data after the rename and
   after the switch.  The sync of a directory is best effort: its
   result is not checked, because not every file system supports it.

A run that failed or was killed leaves <root>/data/.<date>.partial
with the structure files it had fetched and the answers to its HEAD
requests.  Just run snapshot.pl again, the same day or later
(--fresh discards the directory instead).  The rerun
  - fetches the holdings files, the sequence files and the index
    files again, and makes the link tree again;
  - sends no HEAD request again for a file whose listed date is
    still the same;
  - keeps a fetched structure file without a request if its bytes
    are known to be the server's: its length is the server's, and
    its MD5 is the ETag, or it gives the multipart ETag, or the
    ETag is the strong ETag it was fetched with and the file still
    has the MD5 it had then.  A file from a server answer without a
    strong ETag is fetched again;
  - downloads again every file that the failed run downloaded only
    to compare and did not keep.
A dated directory therefore never exists in an incomplete state.

Output is a table of what will be done, the guards, and a summary
(fetched, reused, checked, removed, failed, repaired) with a note for
everything unusual.  --verbose adds a line per file.  Exit status 0
means the snapshot was published.

Removed files.  files/<kind>/removed lists, as before, the names of
files that were in the previous snapshot and are not in this one.
ParsePDB sets the obsolete date of every entry named in
files/hash/removed and files/obs-hash/removed.  Therefore a file
whose entry is still listed (it only lost its PDB-format file, say)
is not written there; it is in files/REMOVED.tsv only, with the
reason "still-listed".  files/REMOVED.tsv has a line for every
removed file: the first ten columns of its manifest line and the
reason (moved: the entry is now in the other holdings file;
still-listed; vanished).  It is how a later run recognises a file
that comes back.

A file that comes back and is the same as the archived copy is
linked to that copy, so nothing appears in files/ of the new
snapshot and ParsePDB does not see it: if the entry's obsolete date
was set when the file was removed, it stays set.  Such files are
listed in files/REAPPEARED.tsv and in the summary, and the run needs
--max-reappeared.  The lasting fix is for the import to take
obsolescence from files/holdings/all_removed_entries.json.gz.

Damaged archive.  A link or file of the previous snapshot that does
not exist is fetched again, up to --max-repaired files; any other
error reading it fails the file.  An archived file that is empty, or
whose size is not the one in the manifest, fails the file.  So does
one whose content is not what the manifest records, but that is only
seen when the file is read: with --verify, or when the run has to
download the file to compare anyway.  An ordinary run does not
notice a damaged file that still has the recorded size.

Run again with --refetch <path> (the path as in the manifest; may be
repeated).  --refetch always downloads the file, whatever state the
archived file is in, and compares the download with the archived
file as it is on disk, not with the manifest:
  - same bytes, or same content compressed differently: the archived
    file is kept, and the new manifest gets its real size and
    checksums (the summary says so when the old manifest was wrong);
  - different, and the archived file does not have the recorded size
    or MD5, or nothing is recorded for it: the download is stored.
    If the manifest has checksums for the file and the download
    agrees with them, the new copy gets the recorded modification
    time, so nothing changes for ParsePDB; otherwise the summary
    says that it may be taken as a new release;
  - different, and the archived file is as recorded: the server's
    file has changed, and the download is stored as changed content.
The damaged file itself is never touched; the new copy is in files/
of the new snapshot.  Older snapshots still link to the damaged
file.

Frozen snapshots.  A snapshot that is the freeze date of a SCOP or
SCOPe release is frozen by hand, by making it read-only
("chmod -R a-w <root>/data/<YYMMDD>"), and is never changed afterwards.
A snapshot counts as frozen when its directory is not writable by its
owner; the check is on the permission bits, so it holds for root too.
  - A run never writes into any earlier snapshot, frozen or not; later
    snapshots only link to the files stored in it.
  - prune.pl and makelinks.pl refuse a frozen snapshot.
  - --manifest-only (with or without --verify) refuses it unless
    --manifest-out names a file outside the snapshot.
There is no option to override this; to change a frozen snapshot it
has to be made writable again by hand.

The manifest.  files/MANIFEST.tsv has one line for each file in the
membership of the snapshot, tab-separated, "-" for unknown:

    path           below <YYMMDD>/snapshot/, e.g. hash/ab/pdb1abc.ent.gz
    entry          PDB id ("-" for sequence and index files)
    kind           pdb, pdb-obs, cif, cif-obs, xml, xml-obs, bundle,
                   seq, index
    storage        the file itself, below <root>/data/
    size           bytes
    md5            of the file as stored
    sha256         of the decompressed content (of the file itself if
                   it is not gzipped)
    last_modified  modification time of the file, UTC; this is the
                   server's Last-Modified when the file was fetched
    etag           last strong ETag seen from the server
    listed_date    date listed for the entry in the holdings
    status         fetched (in this snapshot) or reused
    source         snapshot whose files directory holds the file

Lines starting with "#" are comments; the last line is "#rows" and
the number of rows, so a truncated manifest is recognised.
makelinks.pl remakes the link tree of a pruned snapshot from its
manifest alone ("makelinks.pl <YYMMDD>").

Checksums of reused files are copied from the previous manifest.
When the previous snapshot has none, the reused files have "-" for
md5 and sha256, and the summary says how many.  To fill them in:

    snapshot.pl --manifest-only <YYMMDD>

This reads every file that has no checksum yet (for the first
snapshot made with PAST 2.0, nearly the whole archive membership), so
run it when the machine is otherwise idle.  It can be stopped and
started again: finished rows are kept in MANIFEST.tsv.partial.  It
fails, and names the files, if a file is missing, has the wrong size
or fails the gzip test, or if a file it reads does not have a
checksum that the manifest already records for it.  It also builds a
manifest for an unpruned snapshot of PAST 1.x.  --manifest-out FILE
writes the manifest (and the work file) elsewhere (old snapshots are
read-only); --reference MANIFEST (may be repeated) takes checksums of
the same stored files from another manifest instead of reading the
files.  --hash-reused makes snapshot.pl itself compute the missing
checksums during a run.  Neither reads a file whose checksums are
already recorded.

To check a snapshot against its manifest:

    snapshot.pl --manifest-only <YYMMDD> --verify

This reads every file in the membership of the snapshot, those of
earlier snapshots included, and compares size, MD5 and the SHA-256
of the content with the manifest.  It fails and names every file
that differs, and writes nothing then.  A complete manifest that
verifies is left as it is; missing checksums are filled in as
without --verify.  It can be stopped and started again in the same
way; a file that an interrupted pass has read is not read again (the
output says how many), so remove MANIFEST.tsv.partial to start from
the beginning.  --reference cannot be used with it.

Configuration.  Defaults are at the top of snapshot.pl.  Each can be
given on the command line, and these also in the environment:

    --root DIR         PAST_ROOT         archive root
    --source URL       PAST_SOURCE       base URL of the wwPDB archive
    --tmpdir DIR       PAST_TMPDIR       scratch space (about 20 MB)
    --onhold-url URL   PAST_ONHOLD_URL   on-hold sequences
    --date YYMMDD      PAST_DATE         name of the snapshot (today)

"snapshot.pl --help" lists all options.  prune.pl takes --root or
PAST_ROOT, makelinks.pl PAST_ROOT.  The file TESTING describes the
test harness.

Not archived: entries whose id is not four characters, files that the
holdings list at an unexpected path, and bundles of removed entries.
Each is counted in the output.
