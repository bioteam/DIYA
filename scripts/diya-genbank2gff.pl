#!/usr/bin/perl -w
# $Id: diya-genbank2gff.pl 219 2008-08-26 14:10:34Z briano $;

use Bio::Perl;
use Bio::Tools::GFF;

my ($INFILE, $OUTFILE) = @ARGV;

my $seqi = Bio::SeqIO->new(
	-file	=> $INFILE
);
my $seq = $seqi->next_seq;

my $gff = Bio::Tools::GFF->new(
	-gff_version	=> 3,
	-file			=> ">$OUTFILE",
);

for my $feature ($seq->get_SeqFeatures()) {
	if ($feature->has_tag('score')) {
		($score) = $feature->get_tag_values('score');
		$feature->remove_tag('score');
		$feature->score($score);
	}
	
	$gff->write_feature($feature);
}
