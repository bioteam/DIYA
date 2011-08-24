#!/usr/bin/perl

use strict;
use lib '/site/perl';
use Bio::SeqIO;

#my $in = Bio::SeqIO->new(-file => "2011_08_24_15_52_05-BDRD::rnammer.out.gbk", 
my $in = Bio::SeqIO->new(-file => "singlescore.gbk", 
								 -format => 'genbank');
my $seq = $in->next_seq;
my $out = Bio::SeqIO->new(-file => ">test-seqio4.gbk", 
								  -format => 'genbank');
$out->write_seq($seq);
