#!/usr/bin/perl -w
# $Id: convert2fasta.pl 203 2008-08-13 20:36:59Z briano $

=head1 Description

A script to convert a Genbank file to a fasta file.

=cut

use strict;
use Bio::SeqIO;
use Getopt::Long;

my ($in,$out);

GetOptions( "i=s" => \$in,
			   "o=s" => \$out );

die "$0: need an input file name and an output file name" 
  unless ( $in && $out );

my $converter = Bio::SeqIO->new(-file => ">$out", -format => 'fasta');

my $importer = Bio::SeqIO->new(-file => $in, -format => 'genbank');
my $seq = $importer->next_seq;

$converter->write_seq($seq);


__END__
