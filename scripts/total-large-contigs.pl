#!/usr/bin/perl
# $Id: total-large-contigs.pl 192 2009-01-16 14:37:11Z briano $

use strict;
use Bio::SeqIO;

my $total = 0;
my $totalBig = 0;
my $cutoff = 250;

my $file = shift or die "Need sequence file";

my $in = Bio::SeqIO->new(-file => $file);

while ( my $seq = $in->next_seq ) {
	my $len = $seq->length;
	$total += $len;
	$totalBig += $len if ( $len > $cutoff );
}

print "$file - Total bp: $total, total in contigs greater than $cutoff: $totalBig\n";
