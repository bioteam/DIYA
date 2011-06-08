# $Id: glimmer3ann.pm 290 2008-12-17 19:41:58Z briano $
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

glimmer3ann.pm

=head1 SYNOPSIS

Do not use this module directly. This module is used by diya.pm when it
runs a pipeline.

=head1 DESCRIPTION

This is a variant of glimmer3.pm that uses the script step 'download'
as input. See examples/download-and-annotate.conf.

The gene prediction component of the Do It Yourself Annotator utilizing glimmer
for ab initio detection of potential coding sequences.  

This script has the glimmer3 package in mind, so use with older versions 
of glimmer may require some tampering.  Speficially, use of shell scripts
in glimmer3 to control the individual processes (longorfs, etc) has been 
used here, as opposed to running the individual glimmer programs.

=head1 AUTHOR

Andrew Stewart, andrew.stewart@med.navy.mil
Brian Osborne, briano@bioteam.net

=head1 CONTRIBUTORS

Tim Read, timothy.read@med.navy.mil

=cut

package diya::glimmer3ann;

use strict;
use vars qw(@ISA);
use diya qw($MYSTRAIN);
@ISA = qw(diya);
use Bio::Tools::Glimmer;
use Bio::SeqIO;

sub parse {
	my ($self,$diya) = @_;

	my $PRIMARY_TAG = "gene";
	my $LOCUS_TAG_PREFIX = "orf"; # locus tag prefix
	my $LOCUS_TAG_NUMBER = 0;		# starting number

	my $out = $diya->_outputfile('download');

	print "Parsing $out.predict\n" if $diya->verbose;

	my $parser = Bio::Tools::Glimmer->new(-file   => "$out.predict",
													  -format => 'Glimmer',
													  -detail => "$out.detail" );

	my $gbk = $MYSTRAIN . ".gbk";
	print "Input sequence file is \'$gbk\'\n" if $diya->verbose;
	my $in = Bio::SeqIO->new(-file => $gbk, -format => 'genbank');
	my $seq = $in->next_seq;

	while ( my $feature = $parser->next_feature ) {
		my %tags;
		$tags{locus_tag}	= $LOCUS_TAG_PREFIX . ($LOCUS_TAG_NUMBER += 10);

		my $feat = Bio::SeqFeature::Generic->new(-primary		=> $PRIMARY_TAG,
															  -source_tag	=> 'glimmer3',
															  -start 		=> $feature->start,
															  -end			=> $feature->end,
															  -strand		=> $feature->strand,
															  -tag			=> { %tags },
															 );
		# Add feature to seq object
		$seq->add_SeqFeature($feat);
	}

	# Sort features by location
	my @features = $seq->get_SeqFeatures;
	# sort features by start position
	@features = sort { $a->start <=> $b->start } @features;
	$seq->remove_SeqFeatures;
	$seq->add_SeqFeature(@features);

	# Output
	my $outfile = $out . ".gbk";
	print "Output Genbank file is \'$outfile\'\n" if $diya->verbose;	

	my $seqo = Bio::SeqIO->new(-format	=> 'genbank',
										-file	=> ">$outfile");
	$seqo->write_seq($seq);
}

1;

__END__
