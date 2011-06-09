#!/usr/bin/perl
# $Id: make-fa-and-gbk.pl 198 2009-01-30 16:16:42Z briano $
# Takes a Genbank file, makes a fasta file from it

use strict;
use Bio::SeqIO;

my $gbk = shift or die "No Genbank file given";

my $in = Bio::SeqIO->new( -file => $gbk,
								  -format => 'genbank' );

my ($genome) =  $gbk =~ /(\S+)\.gbk/;

my $out = Bio::SeqIO->new( -file => ">$genome.fa",
								   -format => 'fasta' );

my $seq = $in->next_seq;

$out->write_seq($seq);
