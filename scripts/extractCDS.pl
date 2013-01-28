#!/usr/bin/perl

=head1 NAME

extractCDS.pl

=head1 DESCRIPTION

Extracts the CDS features from a GenBank file and creates nucleotide and
peptide fasta files containing all genes and translations.

     gene            complement(32855..33520)
                     /locus_tag="test0001_380"
                     /inference="ab initio prediction:glimmer3:3.0.2"
     gene            33637..33936
                     /locus_tag="test0001_390"
                     /inference="ab initio prediction:glimmer3:3.0.2"

=cut

use strict;
use Bio::SeqIO;
use Bio::Seq;
use Data::Dumper;

my ($inputfile,$outputfile) = @ARGV;

die "Need both input file name and output file name" 
  unless ( $inputfile && $outputfile );

my $in = Bio::SeqIO->new( -file => "$inputfile.gbk",
			  -format => 'genbank' );

my $ntout = Bio::SeqIO->new( -file => ">>$outputfile",
			   -format => 'fasta' );

my $pepout = Bio::SeqIO->new( -file => ">>$outputfile.pep",
			      -format => 'fasta' );

my $gbk = $in->next_seq;

for my $feature ( $gbk->get_SeqFeatures ) {
    if ( $feature->primary_tag eq "gene" ) {
		my @ids = $feature->get_tag_values("locus_tag");
		my $str = $feature->seq->seq;
		my $seq = Bio::Seq->new( -seq => $str,
								 -display_id => "@ids" );
		$ntout->write_seq($seq);

		$seq = $seq->revcom if ( $feature->location->start > $feature->location->end );
		my $pep = $seq->translate(-complete => 1);
		$pepout->write_seq($pep);
    }
}
