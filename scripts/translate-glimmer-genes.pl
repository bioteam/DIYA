#!/usr/bin/perl -w

use strict;
use Bio::Tools::Glimmer;
use Bio::SeqIO;
use Data::Dumper;

my $MYSEQID = shift or die "No id given";
my $count;

print "Parsing $MYSEQID.predict\n";

my $parser = Bio::Tools::Glimmer->new(-file => "$MYSEQID.predict");
#													 -format => 'Glimmer');
my $pfa = "$MYSEQID-proteins.fa";

my $pout = Bio::SeqIO->new(-file => ">$pfa",
								 -format => 'fasta' );

while ( my $gene = $parser->next_prediction ) {
	$count++;
	#$pout->write_seq($gene);
	#print "Wrote gene " . $gene->display_id . "\n";
}

print "$count predicted genes\n";

__END__
