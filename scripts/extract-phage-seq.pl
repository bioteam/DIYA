#!/arch/bin/perl
=head1 NAME

extract-phage-seq.pl

=head1 DESCRIPTION

Takes a sequence file and phage_finder output, extracts prophage sequence from
the input and writes them to an output fasta file. Usage:

 extract-phage-seq.pl -f IFC.fa -p 2011_09_09_07_32_47-BDRD::phagefinder.out

This command will create a sequence file called IFC_phagefinder.fa.

=cut

use strict;
use Bio::SeqIO;
use Bio::Seq;
use Getopt::Long;

my ( $fasta, $pfoutput, @coords );

GetOptions( "f=s" => \$fasta, "p=s" => \$pfoutput );

die "Need fasta file and phage_finder output file"
  if ( !$fasta || !$pfoutput );

open MYIN, $pfoutput or die "Cannot open file $pfoutput";

my ($id)	= $fasta =~ /^([^.]+)/;
die "Cannot get id from file name $fasta" if ! $id;

# Get phage_finder coordinates and name information
while (<MYIN>) {
    my @line = split /\t/;
    push @coords,
        ( $line[3] - 1 ) . ','
      . ( $line[4] - $line[3] + 1 ) . ','
      . $line[3] . '-'
      . $line[4] . ' '
      . $line[6] . ' '
		. $line[7] . ' '
      . $line[12] . ' '
      . $line[13];
}

my $in = Bio::SeqIO->new(
    -file   => "$fasta",
    -format => 'fasta'
);
my $sequence = $in->next_seq->seq;

my $output = $id . '_phagefinder.fa';
my $out = Bio::SeqIO->new(
    -file   => ">$output",
    -format => 'fasta'
);

my $count = 1;
for my $coord (@coords) {
    my ( $offset, $len, $header ) = split /,/, $coord;
    my $substr = substr( $sequence, $offset, $len );
	 my $name = $id . '_phage_' . $count++; 
    my $seq = Bio::Seq->new( -seq => $substr, -display_id => $name, -desc => $header );
    $out->write_seq($seq);
}
