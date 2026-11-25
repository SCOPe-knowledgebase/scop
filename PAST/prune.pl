#!/lab/bin/perl
#
#   prune.pl - pruning back an old PDB snapshot to save on inodes.
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
# Description
#
# This program will prune back an old snapshot of the pdb to
# save on inodes... the "snapshot" directory is deleted, but the
# "files" will still be there.  All info can be reconstructed
# with "makelinks.pl"
#
# The argument must be the six-digit name of a snapshot (YYMMDD) and
# nothing else: no path, no trailing slash, not "latest".  The snapshot
# that "latest" points to is never pruned.
#
# Needs only modules that come with Perl (Cwd, File::Path, Getopt::Long).
#
# The root of the pdb directory is taken from --root, or the environment
# variable PAST_ROOT, or the default below.

use strict;
use warnings;
use Cwd qw(realpath);
use File::Path qw(remove_tree);
use Getopt::Long;

# ################## config here ########################

# where is pdb directory on this machine?

my $root = $ENV{PAST_ROOT} // "/lab/db/pdb";

# ################## end of config section ###############

# A snapshot whose directory is not writable by its owner is frozen (it is
# the freeze of a SCOP/SCOPe release) and is never changed, whoever runs
# this; -w alone would let root through.
sub frozen {
    my ($d) = @_;
    my @st = stat($d);
    return (! @st) || (! ($st[2] & 0200)) || (! -w $d);
}

MAIN : {
    GetOptions('root=s' => \$root) or die("syntax:  prune.pl [--root DIR] snapshot\n");
    if (@ARGV != 1) {
        die("syntax:  prune.pl [--root DIR] snapshot\n");
    }

    my $name = $ARGV[0];
    if ($name !~ /\A[0-9]{6}\z/) {
        die("'$name' is not a snapshot name (six digits, YYMMDD); nothing pruned.\n");
    }

    my $data = "$root/data";
    my $dir = "$data/$name";

    if (! -d $data) {
        die("root directory $data does not exist!\n");
    }
    if (-l $dir || ! -d $dir) {
        die("snapshot '$name' does not exist (or is not a plain directory)!\n");
    }
    if (! -l "$data/latest") {
        die("'latest' link does not exist!\n");
    }

    my $realLatest = realpath("$data/latest");
    my $realDir = realpath($dir);
    if (! defined $realLatest || ! -d $realLatest) {
        die("'latest' does not point to a directory; nothing pruned.\n");
    }
    if (! defined $realDir) {
        die("cannot resolve '$dir'; nothing pruned.\n");
    }
    if ($realDir eq $realLatest) {
        die("you can't prune the latest snapshot.\n");
    }
    if ($realDir ne realpath($data)."/$name") {
        die("'$name' does not resolve to a directory in $data; nothing pruned.\n");
    }

    my $target = "$dir/snapshot";
    if (-l $target) {
        die("'$target' is a symlink; nothing pruned.\n");
    }
    if (! -d $target) {
        die("snapshot '$name' has no link tree (already pruned?).\n");
    }
    if (! -d "$dir/files") {
        die("snapshot '$name' has no files directory; nothing pruned.\n");
    }
    if (frozen($dir)) {
        die("snapshot '$name' is frozen (read-only); it is never pruned.\n");
    }

    remove_tree($target, { safe => 1, error => \my $err });
    if (@$err || -e $target) {
        die("could not remove all of $target\n");
    }
}
