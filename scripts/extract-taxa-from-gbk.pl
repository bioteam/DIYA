#!/usr/bin/perl -w

=head1 NAME

extract-taxa-from-gbk.pl

=head1 DESCRIPTION

Creates a MEGAN format file from a Genbank file created by Diya

=head1 USAGE

 ./extract-taxa-from-gbk  -f test.gb -i someID -d somedirectory

=cut

use strict;
use lib '/site/perl';
use Getopt::Long;
use Bio::SeqIO;

my ($gbk,$id,$taxa,$outputdir);

GetOptions( "f=s" => \$gbk,
				"i=s" => \$id,
				"d=s" => \$outputdir );

die "No file specified with -f" if ! $gbk;
die "No id specified with -i" if ! $id;

my @cds = grep { $_->primary_tag eq 'CDS' } Bio::SeqIO->new(-file => $gbk)->next_seq->get_SeqFeatures;

# /product="Molybdopterin-guanine dinucleotide biosynthesis
# protein B n=28 Tax=Enterobacteriaceae RepID=A9QYH4_YERPG"

for my $cds ( @cds ) {
	if ( $cds->has_tag('product') ) {
		for my $product ( $cds->get_tag_values('product') ){
			my ($taxon) = $product =~ /Tax=(.+?)RepID/s;
			$taxon =~ s/\s+$//;
			$taxa->{$taxon}++;
		}
	}
}

open MYIN,">$id.csv";

for my $taxon ( keys %{$taxa} ) {
	print MYIN "$taxon," . $taxa->{$taxon} . "\n";
}

`cp $id.csv $outputdir` if ( -d $outputdir );

__END__
