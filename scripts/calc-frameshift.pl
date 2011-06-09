#!/usr/bin/perl 
# $Id: calc-frameshift.pl 170 2008-09-12 16:47:11Z briano $

=head1

Read the "discrp" file, figure out how long the overlap is.

"DiscRep:OVERLAPPING_GENES::975 genes overlap another gene on the same strand."

Example lines:

    Geneyberc_440lcl|contig01102:21475-23070yberc_440
    Geneyberc_450lcl|contig01102:23067-24623yberc_450

=cut

    use strict;

my ($start, $end, $gene);

my $file = shift or die "No discrp file";

open MYIN,$file;

while (<MYIN>) {
    my @line = split,/\s+/;
    $line[2] =~ /(\d+)-(\d+)$/;
    
    if ( $start && $end && $2 && $1 ) {
	my ($fs,$overlap) = get_overlap($start, $end, $1, $2);
	if ( defined $overlap && defined $fs ) {
	    my ($locus1) = $gene =~ /_(\d+)$/;
	    my ($locus2) = $line[1] =~ /_(\d+)$/;
	    
	    # print if the loci are next to each other
	    print "$gene,$line[1]\t$overlap\t$fs\n" 
		if ( $locus1 + 10 == $locus2 );
	}
    }

    $gene = $line[1];
    ($start,$end) = ($1,$2);
}

sub get_overlap {
    my ($start1,$end1,$start2,$end2) = @_;
    my ($end,$begin);

    if ($start1 < $end1) {
	$end = $end1; 
    } else {
	$end = $start1;
    } 

    if ($start2 < $end2) {
	$begin = $start2;
    } else {
	$begin = $end2; 
    } 

    return undef if ( $begin > $end );

    my $fs = ($end - $begin + 1) % 3;
    my $overlap = $end - $begin + 1;
    ($fs, $overlap);
}
