#!/usr/bin/perl -w
# $Id: download-genome.pl 188 2008-08-12 15:22:09Z briano $

=head1 Description

A script to download a fasta file of a genome given an NCBI
project id and an output file name. Example usage:

  download-genome.pl -i 3 -f Borrelia.fa

In this example the script will download the nucleotide sequence for 
NCBI project #3 and write it to the file 'Borrelia.fa'.

=cut

use strict;
use Bio::DB::GenBank;
use Bio::DB::Query::GenBank;
use Getopt::Long;

my ($id,$file);

GetOptions( "i=i" => \$id,
			   "f=s" => \$file );

die "$0: need an NCBI project id and an output file name" 
  unless ( $id && $file );

my $db = Bio::DB::GenBank->new();

my $query = Bio::DB::Query::GenBank->new(
			-db	  => 'nucleotide',
			-query  => $id . ' [genome project] AND srcdb_refseq[prop]' );

my $stream = $db->get_Stream_by_query($query);			

my $seqo = Bio::SeqIO->new( -format => 'fasta', 
									 -file   => ">$file" );

while ( my $seq = $stream->next_seq ) {			
	$seqo->write_seq($seq);
}

__END__
