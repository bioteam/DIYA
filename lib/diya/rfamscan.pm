# $Id: rfamscan.pm 270 2008-12-09 14:02:48Z briano $
#--------------------------------------------------------------------------
# ©Copyright 2008
#
# This file is part of DIYA.
#
# DIYA is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# DIYA is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with the diya package.  If not, see <http://www.gnu.org/licenses/>.
#--------------------------------------------------------------------------

=head1 NAME

rfamscan.pm

=head1 SYNOPSIS

Do not use this module directly. This module is used by diya.pm when it
runs a pipeline.

=head1 DESCRIPTION

A parser for rfam_scan.pl.

The XML could look something like this:

  <parser>
    <executable>rfam_scan.pl</executable>
    <home>/usr/local/bin</home>
    <command>-f gff -o OUTPUTFILE INPUTFILE</command>
    <name>rfam_scan</name>
    <inputformat>fasta</inputformat>
    <inputfrom></inputfrom>
  </parser>


=head1 AUTHOR

Brian Osborne, briano@bioteam.net

=head1 CONTRIBUTORS

Tim Read, timothy.read@med.navy.mil
Andrew Stewart, andrew.stewart@med.navy.mil

=cut

package diya::rfamscan;

use strict;
use base 'diya';
use Bio::Tools::GFF;

sub parse {
	my ($self,$diya) = @_;

	my $rfamout = $diya->_outputfile("rpsblastCDS");

	my $gff = Bio::Tools::GFF->new( -file => $rfamout, 
											  -gff_version => 3 );

	my $seq = $diya->_sequence;

	while (my $line = $gff->next_feature) {
		$line->primary_tag('ncRNA');
		$seq->add_SeqFeature( $line );
	}

	my $seqo = Bio::SeqIO->new(-file	=> ">$rfamout.gbk",
										-format	=> 'genbank');

	$seqo->write_seq($seq);

}


__END__

 GetOptions(	"outfile=s"				=> \$OUTFILE,
					"rfam=s"				=> \$RFAMPATH,
					"bin=s"					=> \$BIN,
			  );

$INFILE = shift @ARGV;
$OUTFILE = $INFILE unless ($OUTFILE);
my $fileroot = $INFILE;
$fileroot =~ s/\.gbk//;

my $seqi = Bio::SeqIO->new( -file => $INFILE );
my $seq = $seqi->next_seq;

unless (-e "$fileroot.fsa") {
	my $tempseq = Bio::SeqIO->new(
		-file	=> ">$fileroot.fsa",
		-format	=> 'fasta');	
	$tempseq->write_seq($seq);
}

my $gff = Bio::Tools::GFF->new( -file => "$fileroot.rfam" , -gff_version => 3);

while (my $line = $gff->next_feature) {
	$line->primary_tag('ncRNA');
	
	$seq->add_SeqFeature( $line );

}

my $seqo = Bio::SeqIO->new(
	-file	=> ">$OUTFILE",
	-format	=> 'genbank');

$seqo->write_seq($seq);


