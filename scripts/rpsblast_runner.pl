#!/usr/bin/perl
#--------------------------------------------------------------------------
# ©Copyright 2008
#
# This file is part of DIYA.
#
# DIYA is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# DIYA is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with the diya package.  If not, see <http://www.gnu.org/licenses/>.
#--------------------------------------------------------------------------

use strict;
use Getopt::Long;
use Bio::SeqIO;
use File::Basename;

$ENV{BLASTMAT} = "/common/data/blastmat" unless defined($ENV{BLASTMAT});

my ($query, $id, $blast_output, $db, $localdata);
my $chunk_size = 1;
Getopt::Long::Configure(("pass_through", "no_auto_abbrev", "no_ignore_case"));
GetOptions("chunk=i"      => \$chunk_size,
           "localdata=s"  => \$localdata,
           "id=i"         => \$id,
           "d=s"          => \$db,
           "i=s"          => \$query,
           "o=s"          => \$blast_output );
my $more_args = join(" ", @ARGV);

# SGE gives us the string "undefined" in the environment variable above if 
# we're not running in a task array.
$id = $ENV{SGE_TASK_ID} if ( ! $id );
$id = 1 if ($id eq "undefined");

my $timestamp = strftime("%Y_%m_%d_%H_%M_%S", localtime);
my $tmp_output= "/tmp/rpsblast-${timestamp}.tmp";
# my $tmp_output= sprintf("/tmp/rpsblast.%05d.tmp", $$);

#
# Fountain of debugging info
#
my $hostname = qx(/bin/hostname);
chomp $hostname;
print STDERR "hostname  = $hostname\n";
print STDERR "query     = $query\n";
print STDERR "chunk     = $chunk_size\n";
print STDERR "id        = $id\n";
print STDERR "more_args = $more_args\n";
print STDERR "temp output = $tmp_output\n";
print STDERR "db to use... = $db\n";

#
# If local data cache is available and specified on the command line,
# use it.
#
if ($localdata) { 
  print STDERR "localdata = $localdata\n"; 
  $db =~ /^(.*)\/([^\/]+)$/;
  my ($blastdb, $target) = ($1, $2);
  print STDERR "Got blastdb = $blastdb, target = $target\n";
  if ((-e $localdata . "/" . $target . ".psq") || 
      (-e $localdata . "/" . $target . ".pal") ||
      (-e $localdata . "/" . $target . ".nal") ||
      (-e $localdata . "/" . $target . ".nsq"))  {
    $db = $localdata . "/" . $target;
    print STDERR "Using $db, since it's local\n";
  }
}

unless (-e $query) {
  die "Unable to locate query file $query\n";
}

#
# Make a private query file with just our inputs.
#
my $queryIO  = new Bio::SeqIO( -format => 'Fasta', -file => $query);
my ($query_local, $query_path, $query_suffix) = fileparse($query);
my $query_fn = sprintf("/tmp/%s.%d", $query_local, $id);
unless (open(QUERY, ">$query_fn")) {
  die "unable to open private query $query_fn for writing\n";
};
my $out = Bio::SeqIO->new(-format => 'Fasta', 
                          -fh     => \*QUERY);
my $counter = 1;
my $seq_count = 1;
while (my $seq = $queryIO->next_seq()) {
  if ($counter == $id) {
    $out->write_seq($seq);
  }
  if ($seq_count % $chunk_size == 0) { $counter++; }
  $seq_count++;
}
close(QUERY);

#
# Construct the BLAST command
#
my $tmp_output_dir;
unless ($tmp_output_dir) { $tmp_output_dir = "rpsblast_tmp"; }
my ($blast_local, $blast_path, $blast_suffix) = fileparse($blast_output);
my $output_fn = sprintf("$tmp_output_dir/%s.%-5d", $blast_local, $id);
my $blast_cmd = "/arch/bin/rpsblast " . 
                "-i $query_fn  " .
                "-o $tmp_output " . 
                "-d $db " .
                $more_args;
print STDERR "BLAST Cmd:  $blast_cmd\n";
system($blast_cmd);

unless (-e $tmp_output_dir) { `mkdir $tmp_output_dir`; }
my $cmd = "/bin/cp $tmp_output $output_fn\n";
print STDERR "$cmd\n";
print STDERR `$cmd`;

#
# Clean up after ourselves.
#
system("rm -f $query_fn");
system("rm -f $tmp_output");

