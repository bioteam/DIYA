#!/usr/bin/perl -w
# $Id: extractCDS.pl 198 2009-01-30 16:16:42Z briano $

=head1 NAME

extractCDS.pl

=head1 DESCRIPTION

Extracts the CDS features from a GenBank file and creates a 
fasta file containing them.

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

my $out = Bio::SeqIO->new( -file => ">>$outputfile",
								   -format => 'fasta' );

my $gbk = $in->next_seq;

for my $feature ( $gbk->get_SeqFeatures ) {
	if ( $feature->primary_tag eq "gene" ) {
		my @ids = $feature->get_tag_values("locus_tag");
		my $str = $feature->seq->seq;

		my $seq = Bio::Seq->new( -seq => $str,
										 -display_id => "@ids" );
		$out->write_seq($seq);
	}
}

__END__

