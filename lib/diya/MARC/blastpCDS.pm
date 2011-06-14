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

blastpCDS.pm

=head1 SYNOPSIS

Do not use this module directly. This module is used by diya.pm when it
runs a pipeline.

=head1 DESCRIPTION

This module uses a Genbank file and blastp output that corresponds
to the CDS regions of the Genbank file and annotates the CDS
features with BLAST match data. The output is an annotated Genbank
file.

=head1 AUTHOR

Brian Osborne, briano@bioteam.net

=head1 CONTRIBUTORS

Tim Read, timothy.read@med.navy.mil
Andrew Stewart, andrew.stewart@med.navy.mil

=cut

package diya::MARC::blastpCDS;

use strict;
use base 'diya';
use Bio::SearchIO;
use Bio::SeqIO;
use Bio::Index::Blast;
use Data::Dumper;

sub parse {
	my ($self,$diya) = @_;

	my $blastout = $diya->_outputfile("MARC::blastpCDS");
	print "Indexing \'" . $blastout . "\'\n" if $diya->verbose;

	my $index = Bio::Index::Blast->new(-filename => "$blastout.index",
				           -write_flag => 1);
	$index->make_index($blastout);

	my $gbk = $diya->_outputfile("MARC::glimmer3");
	my $in = Bio::SeqIO->new(-file => "$gbk.gbk", -format => 'genbank');
	my $seq = $in->next_seq;

	# Remove all features but add them back after
	# adding match data to the CDS features
	my @feats = $seq->remove_SeqFeatures();

	for my $feat (@feats) {

		if ( $feat->primary_tag eq 'gene' ) {

			my @locus = $feat->get_tag_values('locus_tag');

			my %tags = ();
			my $fh;

			eval { $fh = $index->get_stream($locus[0]); };

			if ( $@ ) {
				print "Problem retrieving " . $locus[0] . " from $blastout.index\n";
				next;
			}

			my $blast_report = Bio::SearchIO->new(-noclose => 1,
															  -format  => 'blast',
															  -fh      => $fh);
			my $result = $blast_report->next_result;

			if ( my $hit = $result->next_hit ) {

				my $hsp = $hit->next_hsp;

				( $tags{locus_tag}, $tags{score}, $tags{product} ) = 
				  ($locus[0], $hsp->bits, $hit->description);

				my $CDSfeat = new Bio::SeqFeature::Generic(-primary => 'CDS',
																	 -start	 => $feat->start,
																	 -end		 => $feat->end,
																	 -strand	 => $feat->strand,
																	 -tag		 => { %tags }	);
				$seq->add_SeqFeature($CDSfeat);
				$seq->add_SeqFeature($feat);

			} else {

				$tags{locus_tag} = $locus[0];
				my $CDSfeat = new Bio::SeqFeature::Generic(-primary => 'CDS',
																	 -start	 => $feat->start,
																	 -end		 => $feat->end,
																	 -strand	 => $feat->strand,
																	 -tag     => { %tags }  );
				$seq->add_SeqFeature($CDSfeat);
				$seq->add_SeqFeature($feat);

			}

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

1;

__END__

