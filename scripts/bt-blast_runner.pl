#!/usr/bin/perl
###################################################################
# This software is proprietary to The BioTeam, Inc.
# This document may not be distributed or duplicated, in any part
# or in any manner without the written permission from The BioTeam,
# except for the express purpose for which it was shared.
# Copyright 2009. All Rights Reserved.
###################################################################
#BSUB -L /usr/bin/perl

=head1 NAME

bt-blast_runner.pl

=head1 DESCRIPTION

The execution portion of bt-batchblast. LSF version.

=cut

use strict;
use Getopt::Long;
use Bio::SeqIO;
use File::Basename;

$ENV{BLASTMAT} = "/Jake/apps/ncbi/data/" if ( ! defined( $ENV{BLASTMAT} ) );

my ( $query, $id, $blast_output, $db, $tmp_output, $localdata );
my $chunk_size = 1;
my $bindir     = '/Jake/apps/bin';

Getopt::Long::Configure(
    ( "pass_through", "no_auto_abbrev", "no_ignore_case" ) );
GetOptions(
    "chunk=i"     => \$chunk_size,
    "localdata=s" => \$localdata,
    "id=i"        => \$id,
    "d=s"         => \$db,
    "i=s"         => \$query,
    "o=s"         => \$blast_output
);
my $more_args = join( " ", @ARGV );

# LSB_JOBINDEX == SGE_TASK_ID
$id = $ENV{LSB_JOBINDEX} unless $id;

$tmp_output = sprintf( "/tmp/blastall.%05d.tmp", $$ );

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
print STDERR "out file  = $tmp_output\n";
print STDERR "db to use = $db\n";

#
# If local data cache is available and specified on the command line,
# use it.
#
if ($localdata) {
    print STDERR "localdata = $localdata\n";
    $db =~ /^(.*)\/([^\/]+)$/;
    my ( $blastdb, $target ) = ( $1, $2 );
    print STDERR "Got blastdb = $blastdb, target = $target\n";
    if (   ( -e $localdata . "/" . $target . ".psq" )
        || ( -e $localdata . "/" . $target . ".pal" )
        || ( -e $localdata . "/" . $target . ".nal" )
        || ( -e $localdata . "/" . $target . ".nsq" ) )
    {
        $db = $localdata . "/" . $target;
        print STDERR "Using $db, since it's local\n";
    }
}

die "Unable to locate query file $query" unless ( -e $query );

#
# Make a private query file in /tmp with just our inputs.
#
my $queryIO = Bio::SeqIO->new( -format => 'Fasta', -file => $query );
my ( $query_local, $query_path, $query_suffix ) = fileparse($query);
my $query_fn = sprintf( "/tmp/%s.%d", $query_local, $id );
print STDERR "input  = $query_fn\n";

die "Unable to open private query $query_fn for writing"
  unless ( open( QUERY, ">$query_fn" ) );

my $out = Bio::SeqIO->new(
    -format => 'Fasta',
    -fh     => \*QUERY
);
my $counter   = 1;
my $seq_count = 1;
while ( my $seq = $queryIO->next_seq() ) {
    if ( $counter == $id ) {
        $out->write_seq($seq);
    }
    $counter++ if ( $seq_count % $chunk_size == 0 );
    $seq_count++;
}
close(QUERY);

#
# Construct the BLAST command
#
my $tmp_output_dir;
$tmp_output_dir = "batchblast_tmp" unless ($tmp_output_dir);
my ( $blast_local, $blast_path, $blast_suffix ) = fileparse($blast_output);
my $output_fn = sprintf( "$tmp_output_dir/%s.%-5d", $blast_local, $id );
my $blast_cmd =
  "$bindir/blastall -i $query_fn -o $tmp_output -d $db $more_args";

print STDERR "BLAST Cmd from $0 on $hostname:  $blast_cmd\n";
system($blast_cmd);

`mkdir $tmp_output_dir` unless ( -e $tmp_output_dir );
my $cmd = "/bin/cp $tmp_output $output_fn\n";
print STDERR "$cmd\n";
print STDERR `$cmd`;

#
# Clean up after ourselves.
#
#unlink($query_fn);
#unlink($tmp_output);
