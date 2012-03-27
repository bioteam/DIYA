#!/usr/bin/perl
###################################################################
# This software is proprietary to The BioTeam, Inc.
# This document may not be distributed or duplicated, in any part
# or in any manner without the written permission from The BioTeam,
# except for the express purpose for which it was shared.
# Copyright 2009. All Rights Reserved.
###################################################################

=head1 NAME

bt-batchrpsblast.pl

=head1 DESCRIPTION

bt-batchrpsblast is a high throughput solution for BLAST processing.  
Its syntax is similar to blastall, except that -o and -i are required arguments
(STDIN and STDOUT streams are not yet supported). LSF version.

=cut

use strict;
use File::Spec;
use Bio::SeqIO;
use Getopt::Long;
use POSIX qw(ceil);
use File::Basename;

my ( $query, $jobname, $sync, $project, $arch );
my $chunk          = 1;
my $output         = "stdout.txt";
my $tmp_output_dir = "rpsblast_tmp";
my $bindir         = '/Jake/apps/DIYA/scripts';

Getopt::Long::Configure( "pass_through", "no_auto_abbrev", "no_ignore_case" );
GetOptions(
    "i=s"       => \$query,
    "o=s"       => \$output,
    "tmp_dir=s" => \$tmp_output_dir,
    "chunk=i"   => \$chunk,
    "sync"      => \$sync,
    "jobname=s" => \$jobname,
    "arch=s"    => \$arch,
    "project=s" => \$project
);
my $more_args = join( " ", @ARGV );

$tmp_output_dir = File::Spec->rel2abs($tmp_output_dir);
$output         = File::Spec->rel2abs($output);
die "ERROR:  Please delete temp directory $tmp_output_dir or specify a new one"
  if ( -e $tmp_output_dir );
`mkdir $tmp_output_dir`;

die "Unable to find query file $query" if ( !-e $query );

# Count the number of sequences in the query file.
my $count = `grep '>' $query | wc -l`;
chomp $count;
print STDERR "Counted $count sequences in query $query\n";
my $splits = ceil( $count / $chunk );
if ( $splits > 500 ) {
    $splits = 500;
    if ( ( $count % 500 ) == 0 ) {
        $splits = 500;
        $chunk  = $count / $splits;
    }
    else {
        $chunk  = int( $count / 500 ) + 1;
        $splits = ceil( $count / $chunk );
    }
}

# Submit an array job of $splits tasks, each doing
# $count / $splits query sequences.
print STDERR "Submitting $splits jobs of $chunk queries.\n";

$jobname = "Job$$-blast" unless $jobname;

print STDERR "BLASTMAT = $ENV{BLASTMAT}\n";
my $cmd = "bsub -o batchrpsblast.stdout -e batchrpsblast.stderr ";

if ( $splits > 1 ) {
    $cmd .= "-J \"$jobname" . "[1-$splits]\" ";
}
else {
    $cmd .= "-J $jobname -K ";
}

if ($project) { $cmd .= "-P $project "; }

$cmd .=
    "$bindir/bt-rpsblast_runner.pl "
  . "-i $query "
  . "-o $output "
  . "--chunk $chunk ";
if ($main::LOCALDATA) { $cmd .= "--localdata $main::LOCALDATA "; }
$cmd .= "$more_args";
print STDERR "bsub cmd:  $cmd\n";
system($cmd);

my ( $base_outputname, $outputpath ) = fileparse($output);

## create cleanup script.
my $cleanfn = "cleanup.sh";
open( OUTPUT, ">$cleanfn" );
print OUTPUT <<EOF;
#!/bin/sh
#BSUB -L /bin/sh

echo "waiting for jobs to finish"

cat $tmp_output_dir/$base_outputname.* > $output
rm $tmp_output_dir/$base_outputname.*
rmdir $tmp_output_dir
echo "Done."
EOF
close OUTPUT;

`chmod +x $cleanfn`;

my $finishUp =
"bsub -J \"$jobname.cleanup\" -o batchrpsblast.stdout -e batchrpsblast.stderr -w $jobname ";
$finishUp .= "-K " if $sync;
$finishUp .= "< $cleanfn";
system($finishUp);

`rm $cleanfn`;

__END__
