#!/usr/bin/perl

# Blast Runner
# Chris Dwan (cdwan@bioteam.net)
# August, 2004
# 
# The execution portion of btbatchblast
#
# But really it's rpsblast.  Ha ha. 
#
#$ -S /usr/bin/perl

use strict;
#use lib qw(/RemotePerl/5.8.6);
use lib qw(/site/perl/lib/perl5);
use Getopt::Long;
use Bio::SeqIO;
use File::Basename;

unless (defined($ENV{BLASTMAT})) {
  $ENV{BLASTMAT} = "/common/data/blastmat";
}

#
# Parse Arguments
#
my ($chunk_size, $query, $id, $blast_output, $db, $tmp_output, $localdata);
$chunk_size = 1;
Getopt::Long::Configure(("pass_through", "no_auto_abbrev", "no_ignore_case"));
GetOptions("chunk=i"      => \$chunk_size,
           "localdata=s"  => \$localdata,
           "id=i"         => \$id,
           "d=s"          => \$db,
           "i=s"          => \$query,
           "o=s"          => \$blast_output
          );
my $more_args = join(" ", @ARGV);

# SGE gives us the string "undefined" in the environment variable above if 
# we're not running in a task array.
unless ($id)         { $id = $ENV{SGE_TASK_ID};    }
if ($id eq "undefined") { $id = 1; }

$tmp_output= sprintf("/tmp/blastall.%05d.tmp", $$);

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
unlink($query_fn);
unlink($tmp_output);

