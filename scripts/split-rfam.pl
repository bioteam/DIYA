#!/usr/bin/perl
# Make individual *aln files from Rfam.full

use strict;

open MYIN,"Rfam.full";
my $text = '';

while (my $line = <MYIN>) {

    if ( $line =~ /^# STOCKHOLM/ && $text ) {
	 my ($id) = $text =~ /(RF\d+)/;
	 open MYALN,">aln/$id.aln";
	 print MYALN $text;
	 $text = '';
    }

    $text .= $line;
}

__END__


