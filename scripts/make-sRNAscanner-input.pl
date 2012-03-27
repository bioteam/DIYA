#!/usr/bin/perl
# Makes an input file for sRNAscanner

use strict;
use warnings;
use Bio::SeqIO;
use Bio::FeatureIO;

my $id        = shift or die "No sequence id given";
my $outputdir = shift or die "No outputdir given";

my $out = Bio::FeatureIO->new(
    -file   => ">>$outputdir/$id.ptt",
    -format => 'ptt'
);

my @cds_features =
  grep { $_->primary_tag eq 'CDS' }
  Bio::SeqIO->new( -file => "$outputdir/$id.gbk" )->next_seq->get_SeqFeatures;

for my $feature (@cds_features) {
    $out->write_feature($feature);
}

my $text =
"-----------------------------------------------------------------------------------
Input Options/Parameters used in sRNAscanner-v1
-----------------------------------------------------------------------------------
(1) Number of input Matrices needs to be analyzed:3
(2) Name of the Matrix 1 file:35box_sRNA.matrix
(3) Name of the Matrix 2 file:10box_sRNA.matrix
(4) Name of the Matrix 3 file:terminator.txt.matrix
(5) Which strand needs to be searched for sRNA signals (p=1 or n=2):1
(6) Enter the genome file name (*.fna):${id}.fa
(7) Search the whole genome (g=1) or intergenic (i=2) region:2
(8) Enter the cut-off value 1:2
(9) Enter the cut-off value 2:2
(10) Enter the cut-off value 3:3
(11) Enter the ptt file name:$outputdir/${id}.ptt
(12) Enter the spacer range 1 (between [-35] & [-10] promoter boxes):12-18
(13) Enter the spacer range 2 (sRNA length):40-350
(14) Do you want individual intergenic hits for further analysis, if yes (1) no (2):2
(15) Enter the Unique hit value:200
(16) Enter the minimum Cumulative Sum of Score (CSS):14
----------------------------------------------------------------------------------------
*Note : Please avoid space after colon while specifying parameters
";

open MYOUT, ">$outputdir/Input.data";
print MYOUT $text;

