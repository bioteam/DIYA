#!/usr/bin/perl

=head1 NAME

create-Inputdata.pl

=head1 DESCRIPTION

Script to create the Input.data file for sRNAscanner. Also copies the input sequence
file to the sRNAscanner dir and creates the *ptt file for sRNAscanner.

=cut

use strict;
use Bio::SeqIO;

my $sRNAscannerdir = '/usr/local/share/apps/sRNAscanner';
my $Inputdata      = 'Input.data';

my $id = shift or die "Must provide a sequence id";
my $gbkinfile = shift or die "Must provide a Genbank file";



open MYIN,$idfile or die "Cannot open file $idfile";



__END__
