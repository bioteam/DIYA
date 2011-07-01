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

blastall.pm

=head1 SYNOPSIS

Do not use this module directly. This module is used by diya.pm when it
runs a pipeline.

=head1 AUTHOR

Brian Osborne, briano@bioteam.net

=head1 CONTRIBUTORS

Tim Read, timothy.read@med.navy.mil
Andrew Stewart, andrew.stewart@med.navy.mil

=cut

package diya::MARC::blastall;

use strict;
use vars qw(@ISA);
use diya qw($OUTPUT);
@ISA = qw(diya);

use Bio::SearchIO;
use Bio::SeqIO;
use Data::Dumper;

sub parse {
	my ($self,$diya) = @_;

	my $blastout = $diya->_outputfile("blastall");
	print "Parsing \'" . $blastout . "\'\n" if $diya->verbose;

	my $blastparser = Bio::SearchIO->new(-file => $OUTPUT,
									     -format => 'blast' );
	my $version = $blastparser->version;
	my $program = $blastparser->program;

	my $seq = $diya->_sequence;

	# Remove all features but add them back after
	# adding match data to the CDS features
	my @feats = $seq->remove_SeqFeatures();

	for my $feat (@feats) {
		if ( $feat->primary_tag eq 'CDS' ) {
	 
			my @locus = $feat->get_tag_values('locus_tag');

			my %tags;
			( $tags{locus_tag}, $tags{score}, $tags{product} ) = 
			   get_match_data($blastparser, $locus[0]);
		
			$tags{inference} = "similar to AA sequence:$program:$version";

			my $feat = new Bio::SeqFeature::Generic(-primary => 'CDS',
																 -start	 => $feat->start,
																 -end		 => $feat->end,
																 -strand	 => $feat->strand,
																 -tag		 => { %tags }	);
			$seq->add_SeqFeature($feat);				
	} else {
			$seq->add_SeqFeature($feat);			
	}
}

	my $outfile = $blastout . ".gbk";
 	my $seqo = Bio::SeqIO->new(-file	=> ">$outfile",
 										-format	=> 'genbank');
 	$seqo->write_seq($seq);
	print "Genbank output file is \'$outfile\'\n" if $diya->verbose;
}

# This is not efficient - what's the best way to get a hit by its name?
sub get_match_data {
	my ($parser,$locus) = @_;

	while( my $result = $parser->next_result ) {
		while( my $hit = $result->next_hit ) {
			return ($locus, $hit->bits, $hit->description) 
			  if ( $result->query_name eq $locus );
		}
	}
}

1;

__END__

