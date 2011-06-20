#!/usr/bin/perl

=head1 NAME

gbconvert.pl

=head1 DESCRIPTION

Script to create the tabular and fasta files required by the NCBI
application tbl2asn, then run tbl2asn to create the ASN.1 file 
that NCBI wants.

=cut

use strict;
use diya::MARC::GenbankConvertUtil;
use Getopt::Long;
use Bio::SeqIO;
use FileHandle;
use Bio::Annotation::Collection;
use Bio::Annotation::Comment;

my ($debug, $help, $project, $spacer_start, $spacer_end,
	 $contig_start, $contig_end, $qual, $agp, $taxid);

GetOptions("project|p=i"  => \$project,
			  "help!"        => \$help,
			  "qual|q=s"     => \$qual,
			  "agp!"         => \$agp,
			  "taxid|t=i"    => \$taxid,
			  "debug!"       => \$debug );

usage() if $help;

my $parser = diya::MARC::GenbankConvertUtil->new(-debug => $debug );

my $infile = shift @ARGV or usage('Need a Genbank format file');

my $in = Bio::SeqIO->new(-file => $infile,-format => 'genbank');

my $seq = $in->next_seq;

my $id = $seq->id;

$parser->id($id);
$parser->taxid($taxid) if $taxid;

my $outdir = $parser->outdir($id);

my $definition = $parser->edit_definition($seq->desc);

$seq->desc($definition);

# make output files
my $tblfn = "$outdir/$id.tbl";
my $outfeat = FileHandle->new(">$tblfn") or warn ("$tblfn: $!");

my $fsafn = "$outdir/$id.fsa";
my $outfsa = Bio::SeqIO->new(-file   => ">$fsafn",
									  -format => 'fasta') or warn ("$fsafn: $!");

$parser->make_namemap($infile);

my @oldFeatures = $seq->remove_SeqFeatures;

# fix all the features
FEATURE:
for my $feature ( @oldFeatures ) {

	my $primary_tag = $feature->primary_tag;

	# Each fasta_record describes a contig,
	# the 'cbt' feature is usually right before the 'fasta_record'
	if ( $primary_tag eq 'cbt' ) {
		($spacer_start, $spacer_end) = ($feature->start, $feature->end);
	}

	elsif ( $primary_tag eq 'fasta_record' ) {

		($contig_start, $contig_end) = ($feature->start, $feature->end);
		my $len = $contig_end - $contig_start;

		my @contig_names = $feature->get_tag_values('name');
		my $contig_name = $contig_names[0];

		my @notes = $feature->remove_tag('note') 
		  if ( $feature->has_tag('note') );

		$contig_name = $parser->number_the_duplicate($contig_name)
		  if ( $parser->is_duplicate_name($contig_name) );

		# write to table
		print $outfeat ">Features $contig_name\n";

		if ( my @sfs = $parser->get_feat_from_adj_contig ) {

			for my $sf ( @sfs ) {

				my @locus = $sf->get_tag_values('locus_tag') if 
				  ( $sf->has_tag('locus_tag') );

				my $sf_end = $sf->end;
				my $sf_start = $sf->start;

				# If the feature goes past the 5' end of the current contig (and 
				# we already know that features in this block go past the 3' end) then
				# we must skip this feature, no way to represent it in the *tbl file
				if ( $sf->end <= $contig_end ) {

					# Example: 444  >1 gene
					print $outfeat join("\t", ( ($sf_end - $contig_start + 1) ), 
										  (">1"), $sf->primary_tag ), "\n";

					print "Stored feature: " . $locus[0] . " " .  $sf->primary_tag . 
					  " sf end:$sf_end, contig start:$contig_start\n" if $parser->debug;

					for my $tag ($sf->get_all_tags) {
						for my $value ( $sf->get_tag_values($tag) ) {
							print $outfeat "\t\t\t$tag\t$value\n";
						}
					}

				} elsif ( $sf->primary_tag eq 'rRNA' ) {

					# print $outfeat join("\t", ( ("<$contig_start") ), 
					#					  (">$contig_end"), $sf->primary_tag ), "\n";

               # Example: <1  >111 gene
               print $outfeat join("\t", ( ("<1") ),
                    (">" . $sf_end - $contig_start + 1), $sf->primary_tag ), "\n";


					print "Spanning feature: " . $locus[0] . " " . $sf->primary_tag . 
					  " sf end:$sf_end, contig start:$contig_start\n" if $parser->debug;

					for my $tag ($sf->get_all_tags) {
						for my $value ( $sf->get_tag_values($tag) ) {
							print $outfeat "\t\t\t$tag\t$value\n";
						}
					}

				} else {
					print "Skipped feature: " . $locus[0] . " " . $sf->primary_tag . " start:$sf_start, end:$sf_end, "
					  . "contig start:$contig_start, contig end:$contig_end\n" if $parser->debug;
				}
			}

			$parser->set_feat_from_adj_contig();
		}


		# get coverage statistic
		my $avg;
		for my $note (@notes) {
			($avg) = $note =~ /coverage\s+=\s+([.\d]+)\s+reads/;
			$parser->readsPerBase($len,$avg) if ($avg);
		}

		# write to fasta file
		my $fasta_header = $definition;
		$fasta_header .= " [note=coverage of this contig is $avg" . 'X]';

		my $str = $seq->subseq($contig_start, $contig_end);
		my $featureSeq = Bio::Seq->new(-display_id => $contig_name,
												 -desc       => $fasta_header,
												 -seq        => $str );

		$outfsa->write_seq($featureSeq);
		$outfsa->flush();

		print "fasta_record\t$contig_name\tlength:$len\n" 
		  if $parser->debug;

	}

	# only submitting annotations for 4 primary features
	elsif ( $primary_tag =~ /^(gene|CDS|tRNA|rRNA|repeat_region)$/ ) {

		# skip any feature containing 'N'
		next FEATURE if ( $feature->seq->seq =~ /N/i );

		my @locus = $feature->get_tag_values('locus_tag');

		# skip a feature that starts in spacer ('cbt')
		# next FEATURE if ( $feature->start >= $spacer_start && 
		#						$feature->start <= $spacer_end );
		# skip a feature that ends outside the contig
		# next FEATURE if ( $feature->end > $contig_end );

		# usually a CDS feature and a gene feature are returned here
		my @fixedFeats = $parser->fix_feature($feature);

		$fixedFeats[0] ? $parser->newFeatures(@fixedFeats) : next FEATURE;

		my ($feat_start, $feat_end) = ($feature->start, $feature->end);

	 FIXEDFEATURE:
		for my $feature ( @fixedFeats ) {

			# the feature is entirely contained in the contig
			if ( $feat_start >= $contig_start && $feat_end <= $contig_end ) {

				if ( $feature->strand eq '1' ) {
					print $outfeat join("\t", ($feat_start - $contig_start + 1), 
				  ($feat_end - $contig_start + 1), $feature->primary_tag ), "\n";
				}
				elsif ( $feature->strand eq '-1' ) {
					print $outfeat join("\t", ($feat_end - $contig_start + 1), 
					 ($feat_start - $contig_start + 1), $feature->primary_tag ), "\n";
				}

				print "Feature " . $locus[0] . " is inside contig\n" if $parser->debug;

				for my $tag ($feature->get_all_tags) {
					for my $value ( $feature->get_tag_values($tag) ) {
						print $outfeat "\t\t\t$tag\t$value\n";
					}
				}

				$parser->lastBase($feat_end) if ( $feat_end > $parser->lastBase );
			}
			# 
			elsif (  $feat_start >= $contig_start && $feat_end > $contig_end ) {

				# This feature begins in the next contig, its 3' end is in the spacer,
				# we store this so that we can retrieve it when we're handling the next contig

				if ( $feature->strand eq undef ) {

					print "$locus[0] feature 5' end is in next contig - feature:$feat_end, ".
					"contig:$contig_end, strand=0\n" if $parser->debug;

					$parser->set_feat_from_adj_contig($feature);
					next FIXEDFEATURE;
				}
				# This feature begins in the next contig, its 3' end is in the spacer,
				# we store this so that we can retrieve it when we're handling the next contig
				elsif ( $feature->strand eq '1' && $feat_start > $contig_end ) {

					print "$locus[0] feature 5' end is in next contig - feature:$feat_end, ".
					"contig:$contig_end,strand=1\n" if $parser->debug;

					$parser->set_feat_from_adj_contig($feature);
					next FIXEDFEATURE;
				} 
				# This starts in the contig, strand=1, and ends in the spacer
				# Example: 200	>1575	gene
				elsif ( $feature->strand eq '1' ) {

					print $outfeat join("\t", ($feat_start - $contig_start + 1), 
					(">" . ($contig_end - $contig_start + 1) ), $feature->primary_tag ), "\n";

					print "$locus[0] feature end is past contig - feature:$feat_start-$feat_end, " .
					"contig:$contig_start-$contig_end, strand=1\n" if $parser->debug;

					for my $tag ($feature->get_all_tags) {
						for my $value ( $feature->get_tag_values($tag) ) {
							print $outfeat "\t\t\t$tag\t$value\n";
						}
					}
				}
				# This feature begins in the next contig, its 3' end is in the spacer,
				# we store this so that we can retrieve it when we're handling the next contig
				elsif ( $feature->strand eq '-1' && $feat_start > $contig_end ) {

					print "$locus[0] feature 5' end is in next contig - feature:$feat_end, contig:$contig_end, "
					  . "strand=-1\n" if $parser->debug;

					$parser->set_feat_from_adj_contig($feature);
					next FIXEDFEATURE;
				} 
				# This starts in the contig, strand=-1, and ends in the spacer
				# Example: <444 222 gene
				elsif ( $feature->strand eq '-1' && $feat_start < $contig_end ) {

					print $outfeat join("\t", ("<" . ($contig_end - $contig_start + 1) ), 
				  ( ($feat_start - $contig_start + 1) ), $feature->primary_tag ), "\n";

					print "$locus[0] feature 3' end is in contig - start:$feat_start, end:$feat_end, " .
					  "contig start:$contig_start, contig end:$contig_end, strand=-1\n" if $parser->debug;

					my $startpos =  $parser->get_start_minus($feat_end,$contig_end);
					print $outfeat "\t\t\tcodon_start\t$startpos\n";

					for my $tag ($feature->get_all_tags) {
						for my $value ( $feature->get_tag_values($tag) ) {
							print $outfeat "\t\t\t$tag\t$value\n";
						}
					}
				}
			}
			# This feature begins in the next contig, its 3' end is in the spacer,
			# but there is no strand information, like an rRNA 'gene'
			elsif ( $feat_start > $contig_end ) {

				print "$locus[0] feature 5' end is in next contig - feature:$feat_end, contig:$contig_end,strand=-1\n"
			  if $parser->debug;

				$parser->set_feat_from_adj_contig($feature);
				next FIXEDFEATURE;
			} 
			# One end of the feature is in the contig, the other end
			# precedes the 5' end of the contig
			elsif (  $feat_start < $contig_start && $feat_end <= $contig_end ) {

				print "$locus[0] feature end is before contig - feature:$feat_end, contig:$contig_end\n"
				  if $parser->debug;

				# Example: <1	497	gene
				if ( $feature->strand eq '1' ) {
					print $outfeat join("\t", ("<1"), ($feat_end - $contig_start + 1), 
											  $feature->primary_tag ), "\n";

					my $startpos = $parser->get_start_plus( ($feat_end - $contig_start + 1) );
					print $outfeat "\t\t\tcodon_start\t$startpos\n";
				}
				# Example: 436	>1	gene
				elsif ( $feature->strand eq '-1' ) {

					print $outfeat join("\t", ( ($feat_end - $contig_start + 1) ), 
							 (">1"), $feature->primary_tag ), "\n";
				}

				for my $tag ($feature->get_all_tags) {
					for my $value ( $feature->get_tag_values($tag) ) {
						print $outfeat "\t\t\t$tag\t$value\n";
					}
				}
			}
			# Don't know anything about this feature so skip it
			else {
				print "Skipping feature " . $locus[0] . " with primary tag of $primary_tag\n" 
				  if $parser->debug;
			}
		}
	}
	$outfeat->flush();
}

# calculate overall coverage and add a 'source' feature -
# if the Genbank file is external there's no coverage data
my $comment = $parser->make_top_comment;

if ( $comment ) {
	my $ann = Bio::Annotation::Comment->new;
	$ann->text($comment);
	my $coll = new Bio::Annotation::Collection;
	$coll->add_Annotation('comment',$ann);
	$seq->annotation($coll);
}

# Sort features by location
my @newFeatures = sort { $a->start <=> $b->start } $parser->newFeatures;
$seq->add_SeqFeature(@newFeatures);

$parser->edit_asn_file($seq->desc);

# This first tbl2asn run creates a "discrp" file that has to be parsed
$parser->run_tbl2asn($comment,1);

# Read the discrp file and write a new *.tbl file, fixing discrepancies
$parser->fix_discrp;

# create a quality file for tbl2asn
$parser->create_qual($qual) if $qual;

$parser->create_agp($infile) if $agp;

# Run again
$parser->run_tbl2asn($comment,2);

$parser->cleanup;

sub usage {
	my ($message) = @_;
	$message = "$0 [-h] gbfile" if (! $message );
	print "$message\n";
	exit(-1);
}

__END__

