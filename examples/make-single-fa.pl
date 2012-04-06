#!/usr/bin/perl
# Takes a Genbank file, makes a fasta file from it

use strict;
use Bio::SeqIO;

my $gbk = shift or die "No Genbank file given";

my $dir = shift or die "No directory given";

my $in = Bio::SeqIO->new( -file => "$dir/$gbk",
								  -format => 'genbank' );

my ($genome) =  $gbk =~ /(\S+)\.gbk/;

my $out = Bio::SeqIO->new( -file => ">$genome.fa",
								   -format => 'fasta' );

my $seq = $in->next_seq;

$out->write_seq($seq);
