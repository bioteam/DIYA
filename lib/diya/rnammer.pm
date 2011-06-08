# $Id: rnammer.pm 261 2008-12-03 15:31:54Z briano $
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

rnammer.pm

=head1 SYNOPSIS

Do not use this module directly. This module is used by diya.pm when it
runs a pipeline.

=head1 DESCRIPTION

An rRNA finding component of the diya project which serves
as a result parser for rnammer.

=head1 AUTHOR

Andrew Stewart, andrew.stewart@med.navy.mil
Brian Osborne, briano@bioteam.net

=cut

package diya::rnammer;

use strict;
use base 'diya';
use Bio::Tools::GFF;
use Bio::SeqIO;


sub parse {
	my ($self,$diya) = @_;

	my $LOCUS_TAG_NUMBER = 0;

	my $gbk = $diya->_sequence;

	my $gff = $diya->outputfile('rnammer');

	# Parse rnammer output
	my $gffparser = Bio::Tools::GFF->new( -file => "$gff", 
													  -gff_version => 2);

	while ( my $line = $gffparser->next_feature ) {	
		my %tag;
		$tag{locus_tag} = $gbk->display_id . "_r" . ($LOCUS_TAG_NUMBER += 10);
		$line->set_attributes(-tag => \%tag);
		$line->source_tag('rnammer');
		$gbk->add_SeqFeature( $line );
	}

	# Sort features by location
	my @features = $gbk->get_SeqFeatures;
	# sort features by start position
	@features = sort { $a->start <=> $b->start } @features;

	$gbk->remove_SeqFeatures;
	$gbk->add_SeqFeature(@features);

	my $seqo = Bio::SeqIO->new(-format => 'genbank',
										-file	  => ">$gff.gbk");
	$seqo->write_seq($gbk);
}

__END__
