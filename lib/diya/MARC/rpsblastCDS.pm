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

rpsblastCDS.pm

=head1 SYNOPSIS

Do not use this module directly. This module is used by diya.pm when it
runs a pipeline.

=head1 DESCRIPTION

This module uses a Genbank file and rpsblast output that corresponds
to the CDS regions of the Genbank file and annotates the CDS
features with rpsblast match data. The output is an annotated Genbank
file.

=head1 FILES TO DOWNLOAD

=over 4

=item pub/mmdb/cddid_all.tbl

=item pub/mmdb/cdd/Cdd_LE.tar.gz

=item genomes/CLUSTERS/CURRENT/PRK/PRK_Clusters.bcp

=item genomes/CLUSTERS/CURRENT/PHA/PHA_Clusters.bcp

=back

Then:

 >cat P*Clusters.bcp > Clusters.bcp

Note that the Cdd_LE files are already formatted, you do not run formatrpsbdb.

=head1 AUTHOR

Brian Osborne, briano@bioteam.net

=head1 CONTRIBUTORS

Tim Read, timothy.read@med.navy.mil
Andrew Stewart, andrew.stewart@med.navy.mil

=cut

package diya::MARC::rpsblastCDS;

use strict;

use vars qw(@ISA);
use diya qw($MYSEQID $MYCLUSTERS $MYCDD);
@ISA = qw(diya);

use Bio::SearchIO;
use Bio::SeqIO;
use Bio::Index::Blast;
use Data::Dumper;

sub parse {
	my ($self,$diya) = @_;

	my $clustermap = load_clusters_table();
	my $cddmap = load_cdd_table();

	print "Cluster map from \'$MYCLUSTERS\' loaded\n" 
	  if ( defined $clustermap && $diya->verbose );
   print "CDD map from \'$MYCDD\' loaded\n" 
     if ( defined $cddmap && $diya->verbose );

	my $blastout = $diya->_outputfile("MARC::rpsblastCDS");
	print "Indexing \'" . $blastout . "\'\n" if $diya->verbose;

	my $index = Bio::Index::Blast->new(-filename => "$blastout.index",
												  -write_flag => 1);
	$index->make_index($blastout);

	my $gbk = $diya->_outputfile("MARC::blastpCDS");
	my $in = Bio::SeqIO->new(-file => "$gbk.gbk", -format => 'genbank');
	my $seq = $in->next_seq;

	# Remove all features but add them back after
	# adding match data to the CDS features
	my @feats = $seq->remove_SeqFeatures();

	for my $feat (@feats) {

		if ( $feat->primary_tag eq 'CDS' ) {

			my @locus = $feat->get_tag_values('locus_tag');
			
			my $fh;

			eval { $fh = $index->get_stream($locus[0]); };

			if ( $@ ) {
				print "Problem retrieving " . $locus[0] . " from $blastout.index\n";
				next;
			}

			my $blast_report = Bio::SearchIO->new(-noclose => 1,
															  -format  => 'blast',
															  -fh      => $fh );
			my $result = $blast_report->next_result;

			if ( $result->num_hits > 0 ) {

				my $tags = get_feat_tags($feat);

				my $hit = get_best_hit($result);
				my $hsp = $hit->next_hsp;

				print "Found rpsblast hit " . $hit->description . " for $locus[0]\n" if $diya->verbose;

				# get domain id from hit
				my ($cid) = $hit->description =~ /^(\w+)/;

				if ( defined $clustermap->{$cid} ) {

					$tags->{'score'} = $hit->significance;
					$tags->{'rps_gi'} = $cid;
					$tags->{'inference'} = 'rpsblast';

					# set cluster
					$tags->{'cluster'} = $clustermap->{$cid}->{'entry'};
					$tags->{'product'} = $clustermap->{$cid}->{'definition'};
					$tags->{'group'}	 = $clustermap->{$cid}->{'group'};

			   } elsif ( defined $cddmap->{$cid} ) {

					$tags->{'cluster'} = $cddmap->{$cid}->{'entry'};
               $tags->{'product'} = $cddmap->{$cid}->{'definition'};

				}

				my @products = $feat->get_tag_values('product') if 
				  ( $feat->has_tag('product') );

				# If this feature is already annotated with a UniRef match
				if ( $products[0] =~ /RepID=/ && defined $tags->{'product'} ) {

					my $note = "similar to " . $tags->{'product'};
					$note .= " " . $tags->{'cluster'};
					$feat->add_tag_value('note',$note);

					$seq->add_SeqFeature($feat);

				} else {

					my $CDSfeat = new Bio::SeqFeature::Generic(-primary => 'CDS',
																			 -start	 => $feat->start,
																			 -end		 => $feat->end,
																			 -strand	 => $feat->strand,
																			 -tag		 => { %{$tags} }	);
					# add xrefs
					if ( $cid && $clustermap->{$cid}->{'xref'} ) {
						for ( @{$clustermap->{$cid}->{'xref'}} ) {
							$CDSfeat->add_tag_value( 'xref' => $_ );
						}
					}

					$seq->add_SeqFeature($CDSfeat);
				}

			} else {

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

sub get_feat_tags {
	my $feat = shift;
	my %tags;

   for my $tag ($feat->get_all_tags) {
		my @vals = $feat->get_tag_values($tag);
		$tags{$tag} = $vals[0];
	}
	\%tags;
}

sub get_best_hit {
	my $result = shift;
	my $best;
	my $bits = 0;

	while ( my $hit = $result->next_hit ) {
		my $hsp = $hit->next_hsp;

		if ( $hsp->bits > $bits ) {
			$bits = $hsp->bits;
			$best = $hit;
		}

	}
	return $best;
}

# parse contents of Clusters.bcp into %cluster-map
sub load_clusters_table {

	local $/ = "////\n";

	my %clustermap;

	open MYIN,$MYCLUSTERS or die "Cannot open file $MYCLUSTERS";

	while (<MYIN>) {

		my $rec = $_;

		my ($top, $bottom) = split(/PROTEINS/, $rec);

		my ($entry) = $top =~ /ENTRY\s+(\w*)/;
		my ($definition) = $top =~ /DEFINITION\s+(.+)/;
		my ($group) = $top =~ /COG_GROUP\s+(.+)/;
		my @xrefs = $top =~ /XREF\s+(\w+)/g;

		my %cluster = ('entry'			=> $entry,
							'definition'	=> $definition,
							'group'			=> $group,
							'xref'			=> \@xrefs,
						  );

		$clustermap{$entry} = \%cluster if ( $entry );

		for my $xref ( @xrefs ) {
			$clustermap{$xref} = \%cluster;
		}
	}
	\%clustermap;
}

# parse contents of cddid_all.tbl file into %cddmap
sub load_cdd_table {

	my %cddmap;

	open MYIN,$MYCDD or die "Cannot open file $MYCDD";

	while ( my $line = <MYIN> ) {
		# Examples:
		# 33678	COG3889	COG3889	Predicted solute binding protein [General function prediction only]	872
		# 33328  COG3525 Chb N-acetyl-beta-hexosaminidase [Carbohydrate transport and metabolism] 732
		if ( $line =~ /^\d+\s+(COG\d+)\s+\S+\s+([^[]+)\s/ ) {
			my ($id,$def) = ($1,$2);
			$def =~ s/\s+$//;
			$def =~ s/\s+/ /g;
			my  %cluster = ('entry'      => $id,
								 'definition' => $def );
			$cddmap{$id} = \%cluster;
		}

		if ( $line =~ /^\d+\s+(CHL\d+)\s+([^;]+?)\s+\d+$/ ) {
         my ($id,$def) = ($1,$2);
         $def =~ s/\s+$//;
         $def =~ s/\s+/ /g;
			my  %cluster = ('entry'      => $id,
								 'definition' => $def );
			$cddmap{$id} = \%cluster;
		}

		if ( $line =~ /^\d+\s+(KOG\d+)\s+KOG\d+\s+KOG\d+,\s+([^,[]+)/ ) {
			my ($id,$def) = ($1,$2);
         $def =~ s/\s+$//;
         $def =~ s/\s+/ /g;
			my  %cluster = ('entry'      => $id,
								 'definition' => $def );
			$cddmap{$id} = \%cluster;
		}

		# 103623	PRK08948	PRK08948	2-octaprenyl-6-methoxyphenyl hydroxylase; Validated	392
		if ( $line =~ /^\d+\s+(PRK\d+)\s+\w+\s+([^;]+)/ ) {
			my ($id,$def) = ($1,$2);
         $def =~ s/\s+$//;
         $def =~ s/\s+/ /g;
			my  %cluster = ('entry'      => $id,
								 'definition' => $def );
			$cddmap{$id} = \%cluster;
		}
	}
	\%cddmap;
}

1;

__END__
