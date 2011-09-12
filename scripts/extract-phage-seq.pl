#!/usr/bin/perl
# Takes a sequence file and phage_finder output, extracts prophage sequence from
# the input and writes them to an output fasta file

use strict;
use Bio::SeqIO;
use Getopt::Long;

my ( $fasta, $pfoutput, $output, @coords );

GetOptions( "f=s" => \$fasta, "p=s" => \$pfoutput, "o=s" => \$output );

die "Need fasta file, phage_finder output file name, and output file name"
  if ( !$fasta || !$pfoutput || !$output );

open MYIN, $pfoutput or die "Cannot open file $pfoutput";

# Get phage_finder coordinates and name information
while (<MYIN>) {
    my @line = split /\t/;
    push @coords,
        ( $line[3] - 1 ) . ','
      . ( $line[4] - $line[3] + 1 ) . ','
      . $line[3] . '-'
      . $line[4] . ' '
      . $line[5] . ' '
      . $line[10] . ' '
      . $line[11];
}

my $in = Bio::SeqIO->new(
    -file   => "$fasta",
    -format => 'fasta'
);
my $sequence = $in->next_seq->seq;

my $out = Bio::SeqIO->new(
    -file   => ">$output",
    -format => 'fasta'
);

for my $coord (@coords) {
    my ( $offset, $len, $name ) = split /,/, $coord;
    my $substr = substr( $sequence, $offset, $len );
    my $seq = Bio::Seq->new( -seq => $substr, -display_id => "$output $name" );
    $out->write_seq($seq);
}

