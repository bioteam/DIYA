#--------------------------------------------------------------------------
# ©Copyright 2011
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

GenbankConvertUtil.pm

=head1 SYNOPSIS

Utility methods for Genbank submissions.

=head1 DESCRIPTION

Code to convert Genbank format files to the format required by the tbl2asn,
the application that creates ASN.1 for Genbank. Also does any QC and text
editing on the CDS and gene features as required by NCBI.

This script requires the tbl2asn program from NCBI.

=head1 AUTHOR

Brian Osborne, briano@bioteam.net

=head1 CONTRIBUTORS

Tim Read, timothy.read@med.navy.mil
Andrew Stewart, andrew.stewart@med.navy.mil
Mike Cariaso, cariaso@bioteam.net

=cut

package diya::MARC::GenbankConvertUtil;

use strict;
use Data::Dumper qw(Dumper);
use Bio::SeqIO;
use FileHandle;
use Date::Format qw(time2str);

my ($spacer_start,    $spacer_end,       $contig_start,
    $contig_end);

sub new {
	my ($module,@args) = @_;
	my $class = ref($module) || $module;
	my $self = {};

	bless($self,$class);

	$self->_initialize(@args);

	$self;
}

sub fixAndPrint {
    my ( $self, $ref, $outfeat, $outfsa, $seq, $definition ) = @_;
    my @oldFeatures = @{$ref};

    # Fix all the features
  FEATURE:
    for my $feature (@oldFeatures) {

        my $primary_tag = $feature->primary_tag;

        # Each fasta_record describes a contig, and
        # the 'cbt' feature is usually right before the 'fasta_record'
        if ( $primary_tag eq 'cbt' ) {
            ( $spacer_start, $spacer_end ) = ( $feature->start, $feature->end );
        }
        elsif ( $primary_tag eq 'fasta_record' ) {

            ( $contig_start, $contig_end ) = ( $feature->start, $feature->end );
            my $len = $contig_end - $contig_start;

            my @contig_names = $feature->get_tag_values('name');
            my $contig_name  = $contig_names[0];

            my @notes = $feature->remove_tag('note')
              if ( $feature->has_tag('note') );

            $contig_name = $self->number_the_duplicate($contig_name)
              if ( $self->is_duplicate_name($contig_name) );

            # Write to the *tbl file
            print $outfeat ">Features $contig_name\n";

            if ( my @sfs = $self->get_feat_from_adj_contig ) {

                for my $sf (@sfs) {

                    my @locus = $sf->get_tag_values('locus_tag')
                      if ( $sf->has_tag('locus_tag') );

                    my $sf_end   = $sf->end;
                    my $sf_start = $sf->start;

          # If the feature goes past the 5' end of the current contig (and
          # we already know that features in this block go past the 3' end) then
          # we must skip this feature, no way to represent it in the *tbl file
                    if ( $sf->end <= $contig_end ) {

                        # Example: 444  >1 gene
                        print $outfeat join( "\t",
                            ( ( $sf_end - $contig_start + 1 ) ),
                            (">1"), $sf->primary_tag ),
                          "\n";

                        print "Stored feature: "
                          . $locus[0] . " "
                          . $sf->primary_tag
                          . " sf end:$sf_end, contig start:$contig_start\n"
                          if $self->debug;

                        for my $tag ( $sf->get_all_tags ) {
                            for my $value ( $sf->get_tag_values($tag) ) {
                                print $outfeat "\t\t\t$tag\t$value\n";
                            }
                        }

                    }
                    elsif ( $sf->primary_tag eq 'rRNA' ) {

               # print $outfeat join("\t", ( ("<$contig_start") ),
               #                     (">$contig_end"), $sf->primary_tag ), "\n";

                        # Example: <1  >111 gene
                        print $outfeat join( "\t",
                            ( ("<1") ),
                            ( ">" . $sf_end - $contig_start + 1 ),
                            $sf->primary_tag ),
                          "\n";

                        print "Spanning feature: "
                          . $locus[0] . " "
                          . $sf->primary_tag
                          . " sf end:$sf_end, contig start:$contig_start\n"
                          if $self->debug;

                        for my $tag ( $sf->get_all_tags ) {
                            for my $value ( $sf->get_tag_values($tag) ) {
                                print $outfeat "\t\t\t$tag\t$value\n";
                            }
                        }

                    }
                    else {
                        print "Skipped feature: "
                          . $locus[0] . " "
                          . $sf->primary_tag
                          . " start:$sf_start, end:$sf_end, "
                          . "contig start:$contig_start, contig end:$contig_end\n"
                          if $self->debug;
                    }
                }

                $self->set_feat_from_adj_contig();
            }

            # Get coverage statistic
            my $avg;
            for my $note (@notes) {
                ($avg) = $note =~ /coverage\s+=\s+([.\d]+)\s+reads/;
                $self->readsPerBase( $len, $avg ) if ($avg);
            }

            my $fasta_header = $definition;

            # Write to fasta file if there's coverage data
            # $fasta_header .= " [note=coverage of this contig is ${avg}X"
            #  if $avg;

            my $str = $seq->subseq( $contig_start, $contig_end );
            my $featureSeq = Bio::Seq->new(
                -display_id => $contig_name,
                -desc       => $fasta_header,
                -seq        => $str
            );

            $outfsa->write_seq($featureSeq);
            $outfsa->flush();

            print "fasta_record\t$contig_name\tlength:$len\n"
              if $self->debug;
        }

        # Only submitting annotations for these primary features
        elsif ( $primary_tag =~ /^(gene|CDS|tRNA|rRNA|repeat_region|ncRNA)$/ ) {

            # Skip any feature with sequence containing 'N'
            next FEATURE if ( $feature->seq->seq =~ /N/i );

            my @locus = $feature->get_tag_values('locus_tag')
              if $feature->has_tag('locus_tag');

            # skip a feature that starts in spacer ('cbt')
            # next FEATURE if ( $feature->start >= $spacer_start &&
            #                       $feature->start <= $spacer_end );
            # skip a feature that ends outside the contig
            # next FEATURE if ( $feature->end > $contig_end );

            # usually a CDS feature and a gene feature are returned here
            my @fixedFeats = $self->fix_feature($feature);

            $fixedFeats[0] ? $self->newFeatures(@fixedFeats) : next FEATURE;

            my ( $feat_start, $feat_end ) = ( $feature->start, $feature->end );

          FIXEDFEATURE:
            for my $feature (@fixedFeats) {

                # the feature is entirely contained in the contig
                if ( $feat_start >= $contig_start && $feat_end <= $contig_end )
                {

                    if ( $feature->strand eq '1' ) {
                        print $outfeat join( "\t",
                            ( $feat_start - $contig_start + 1 ),
                            ( $feat_end - $contig_start + 1 ),
                            $feature->primary_tag ),
                          "\n";
                    }
                    elsif ( $feature->strand eq '-1' ) {
                        print $outfeat join( "\t",
                            ( $feat_end - $contig_start + 1 ),
                            ( $feat_start - $contig_start + 1 ),
                            $feature->primary_tag ),
                          "\n";
                    }

                    print "Feature " . $locus[0] . " is inside contig\n"
                      if $self->debug;

                    for my $tag ( $feature->get_all_tags ) {
                        for my $value ( $feature->get_tag_values($tag) ) {
                            print $outfeat "\t\t\t$tag\t$value\n";
                        }
                    }

                    $self->lastBase($feat_end)
                      if ( $feat_end > $self->lastBase );
                }

                elsif ($feat_start >= $contig_start
                    && $feat_end > $contig_end )
                {

  # This feature begins in the next contig, its 3' end is in the spacer,
  # we store this so that we can retrieve it when we're handling the next contig

                    if ( $feature->strand eq undef ) {

                        print
"$locus[0] feature 5' end is in next contig - feature:$feat_end, "
                          . "contig:$contig_end, strand=0\n"
                          if $self->debug;

                        $self->set_feat_from_adj_contig($feature);
                        next FIXEDFEATURE;
                    }

  # This feature begins in the next contig, its 3' end is in the spacer,
  # we store this so that we can retrieve it when we're handling the next contig
                    elsif ($feature->strand eq '1'
                        && $feat_start > $contig_end )
                    {

                        print
"$locus[0] feature 5' end is in next contig - feature:$feat_end, "
                          . "contig:$contig_end,strand=1\n"
                          if $self->debug;

                        $self->set_feat_from_adj_contig($feature);
                        next FIXEDFEATURE;
                    }

                   # This starts in the contig, strand=1, and ends in the spacer
                   # Example: 200  >1575   gene
                    elsif ( $feature->strand eq '1' ) {

                        print $outfeat join( "\t",
                            ( $feat_start - $contig_start + 1 ),
                            ( ">" . ( $contig_end - $contig_start + 1 ) ),
                            $feature->primary_tag ),
                          "\n";

                        print
"$locus[0] feature end is past contig - feature:$feat_start-$feat_end, "
                          . "contig:$contig_start-$contig_end, strand=1\n"
                          if $self->debug;

                        for my $tag ( $feature->get_all_tags ) {
                            for my $value ( $feature->get_tag_values($tag) ) {
                                print $outfeat "\t\t\t$tag\t$value\n";
                            }
                        }
                    }

  # This feature begins in the next contig, its 3' end is in the spacer,
  # we store this so that we can retrieve it when we're handling the next contig
                    elsif ($feature->strand eq '-1'
                        && $feat_start > $contig_end )
                    {

                        print
"$locus[0] feature 5' end is in next contig - feature:$feat_end, contig:$contig_end, "
                          . "strand=-1\n"
                          if $self->debug;

                        $self->set_feat_from_adj_contig($feature);
                        next FIXEDFEATURE;
                    }

                  # This starts in the contig, strand=-1, and ends in the spacer
                  # Example: <444 222 gene
                    elsif ($feature->strand eq '-1'
                        && $feat_start < $contig_end )
                    {

                        print $outfeat join( "\t",
                            ( "<" . ( $contig_end - $contig_start + 1 ) ),
                            ( ( $feat_start - $contig_start + 1 ) ),
                            $feature->primary_tag ),
                          "\n";

                        print
"$locus[0] feature 3' end is in contig - start:$feat_start, end:$feat_end, "
                          . "contig start:$contig_start, contig end:$contig_end, strand=-1\n"
                          if $self->debug;

                        my $startpos =
                          $self->get_start_minus( $feat_end, $contig_end );
                        print $outfeat "\t\t\tcodon_start\t$startpos\n";

                        for my $tag ( $feature->get_all_tags ) {
                            for my $value ( $feature->get_tag_values($tag) ) {
                                print $outfeat "\t\t\t$tag\t$value\n";
                            }
                        }
                    }
                }

          # This feature begins in the next contig, its 3' end is in the spacer,
          # but there is no strand information, like an rRNA 'gene'
                elsif ( $feat_start > $contig_end ) {

                    print
"$locus[0] feature 5' end is in next contig - feature:$feat_end, contig:$contig_end,strand=-1\n"
                      if $self->debug;

                    $self->set_feat_from_adj_contig($feature);
                    next FIXEDFEATURE;
                }

                # One end of the feature is in the contig, the other end
                # precedes the 5' end of the contig
                elsif ($feat_start < $contig_start
                    && $feat_end <= $contig_end )
                {

                    print
"$locus[0] feature end is before contig - feature:$feat_end, contig:$contig_end\n"
                      if $self->debug;

                    # Example: <1   497 gene
                    if ( $feature->strand eq '1' ) {
                        print $outfeat join( "\t",
                            ("<1"), ( $feat_end - $contig_start + 1 ),
                            $feature->primary_tag ),
                          "\n";

                        my $startpos = $self->get_start_plus(
                            ( $feat_end - $contig_start + 1 ) );
                        print $outfeat "\t\t\tcodon_start\t$startpos\n";
                    }

                    # Example: 436  >1  gene
                    elsif ( $feature->strand eq '-1' ) {

                        print $outfeat join( "\t",
                            ( ( $feat_end - $contig_start + 1 ) ),
                            (">1"), $feature->primary_tag ),
                          "\n";
                    }

                    for my $tag ( $feature->get_all_tags ) {
                        for my $value ( $feature->get_tag_values($tag) ) {
                            print $outfeat "\t\t\t$tag\t$value\n";
                        }
                    }
                }

                # Don't know anything about this feature so skip it
                else {
                    print "Skipping feature "
                      . $locus[0]
                      . " with primary tag of $primary_tag\n"
                      if $self->debug;
                }
            }
        }
        $outfeat->flush();
    }
}

sub list_primary_tags {
	my $self = shift;
	my $file = $self->{file};

	my $in = Bio::SeqIO->new(-file => $file,
						     -format => 'genbank');
	my $seq = $in->next_seq;

	my %tags;

	for my $feat ($seq->get_SeqFeatures) {
		my $tag = $feat->primary_tag;
		$tags{$tag}++;
   }

	for my $tag (keys %tags) {
		print "$tag\t$tags{$tag}\n";
	}
}

sub fix_feature {
    my ( $self, $feat ) = @_;

    if ( $feat->primary_tag eq 'CDS' ) {

        # Edit product tag and note tag
        if ( $feat->has_tag('product') ) {

            my ($product) = $feat->remove_tag("product");

            # If 'note' and 'product' are duplicated
            if ( $feat->has_tag('note') ) {
                my @notes = $feat->remove_tag('note');
                for my $note (@notes) {
                    $note =~ s/\s+$//;
                    my $str = "similar to " . $product;
                    $feat->add_tag_value( 'note', $note )
                      unless ( $note eq $str );
                }
            }

            print "Product name before correction: $product\n" if $self->debug;

            $product = trim($product);

            ( $product, $feat ) = fix_uniref( $product, $feat );

            ( $product, $feat ) = fix_cog( $product, $feat );

            ( $product, $feat ) = remove_species( $product, $feat );

            ( $product, $feat ) = move_species( $product, $feat );

            $product = correct_spelling($product);

            $product = add_trailing($product);

            $product = remove_trailing($product);

            $product = make_singular($product);

            $product = remove_banned($product);

            if ( is_hypothetical($product) ) {
                $product = 'hypothetical protein';
                $feat->add_tag_value( 'note', $product );
            }

            $product = remove_similar($product);

            $product = remove_loci($product);

            $product = remove_trailing($product);

            ( $product, $feat ) = remove_semicolon( $product, $feat );

            $feat = edit_note($feat);

            # Finally add it back
            $feat->add_tag_value( 'product', $product );

            print "Product name after correction: $product\n" if $self->debug;
        }
        else {
            # If there is no 'product' then it's a 'hypothetical protein'
            $feat->add_tag_value( 'product', 'hypothetical protein' );
        }

        # Change locus_tag to protein_id
        if ( $feat->has_tag('locus_tag') ) {
            my ($locus) = $feat->remove_tag('locus_tag');
            my $protein_id = 'gnl|' . $self->accession_prefix . '|' . $locus;
            $feat->add_tag_value( 'protein_id', $protein_id );
        }

        # Change xref tags to note tags
        if ( $feat->has_tag('xref') ) {
            my @xrefs = $feat->remove_tag('xref');

            for my $xref (@xrefs) {
                my $note = "similar to $xref";
                $feat->add_tag_value( 'note', $note ) if ( !$xref eq 'family' );
            }
        }

        $feat = fix_rps($feat);

        $feat = remove_scores($feat);

        return $feat;
    }

    elsif ( $feat->primary_tag eq 'gene' ) {

        return $feat;
    }

    elsif ( $feat->primary_tag eq 'repeat_region' ) {

        return $feat;
    }

    elsif ( $feat->primary_tag eq 'ncRNA' ) {

        return $feat;
    }

    elsif ( $feat->primary_tag eq 'tRNA' ) {

        # Add inference tag
        $feat->add_tag_value( 'inference', 'profile:tRNAscan-SE:1.23' )
          if ( !$feat->has_tag("inference") );

        # Edit note tag
        if ( $feat->has_tag("note") ) {
            my ($note) = $feat->remove_tag("note");
            $note =~ s/score=/tRNAscan-SE score=/;
            $feat->add_tag_value( 'note', $note );
        }

        # Change Codon tag to a note tag
        if ( $feat->has_tag("Codon") ) {
            my @codon = $feat->remove_tag("Codon");
            $feat->add_tag_value( 'note', "Anticodon is $codon[0]" );
        }

        # Remove AminoAcid and Name tags
        $feat->remove_tag("AminoAcid") if $feat->has_tag('AminoAcid');
        $feat->remove_tag("Name")      if $feat->has_tag('Name');

        # Remove 'ID', which will be added back as 'product'
        if ( $feat->has_tag("ID") ) {
            my ($ID) = $feat->remove_tag("ID");
            $ID =~ s/:/-/;
            $feat->add_tag_value( 'product', $ID );
        }

        # Add a 'gene' feature for the tRNA
        my $genefeat = Bio::SeqFeature::Generic->new( -primary_tag => 'gene' );
        $genefeat->location( $feat->location );
        $genefeat->strand( $feat->strand );

        my ($loci) = $feat->remove_tag('locus_tag')
          if $feat->has_tag('locus_tag');
        $genefeat->add_tag_value( 'locus_tag', $loci );

        $feat = remove_scores($feat);

        return ( $feat, $genefeat );
    }

    elsif ( $feat->primary_tag eq 'rRNA' ) {

        # Delete note tag with "frame=." value
        if ( $feat->has_tag('note') ) {
            my @notes = $feat->remove_tag("note");

            for my $note (@notes) {
                $note =~ s/score=/RNAMMER score=/;
                $feat->add_tag_value( 'note', $note )
                  if ( $note !~ /frame=\./ );
            }
        }

        # Add a 'gene' feature for the rRNA
        my $genefeat = Bio::SeqFeature::Generic->new( -primary_tag => 'gene' );
        $genefeat->location( $feat->location );
        $genefeat->strand( $feat->strand );

        my ($loci) = $feat->remove_tag('locus_tag')
          if $feat->has_tag('locus_tag');

        $genefeat->add_tag_value( 'locus_tag', $loci );

        $feat = remove_scores($feat);

        return ( $feat, $genefeat );
    }

}

sub remove_scores {
    my $feat = shift;

    # Remove all score tags
    if ( $feat->has_tag('score') || $feat->has_tag('Score') ) {
        my @scores = $feat->remove_tag('score');
    }

    # Remove all score values
    for my $tag ( $feat->get_all_tags ) {
        my @values = $feat->get_tag_values($tag);
        $feat->remove_tag($tag);
        for my $value (@values) {
            $feat->add_tag_value( $tag, $value )
              if ( $value !~ /^\s*score\s*[=:]/i );
        }
    }
    $feat;
}

sub move_species {
    my ( $product, $feat ) = @_;

    # If the product looks something like:
    # "ABC protein [Yersinia pestis CO92]." then correct
    if ( $product =~ /(.+)\s+(\[.+?\])\.$/ ) {

        my ( $newProduct, $species ) = ( $1, $2 );

        # Add back the product but do not use locus tags from other
        # genomes, e.g. product="hypothetical protein YPO0973
        if ( $newProduct =~ /^hypothetical protein\s+.*/i ) {
            $product = 'hypothetical protein';
        }
        else {
            $product = $newProduct;
        }

        my @score = $feat->remove_tag('score') if $feat->has_tag('score');

        if ( $score[0] ) {
            my $note = "similar to $newProduct $species";
            $feat->add_tag_value( 'note', $note );

            #$note = "BLAST score=" . $score[0];
            #$feat->add_tag_value('note',$note);
        }
    }

    ( $product, $feat );
}

sub fix_rps {
    my $feat = shift;

    # Features with 'rps_gi' tag need to be rearranged
    if ( $feat->has_tag('rps_gi') ) {
        my @ids = $feat->remove_tag('rps_gi');

        for my $id (@ids) {
            my $note = "similar to motif $id";
            $feat->add_tag_value( 'note', $note ) unless ( $id eq 'family' );
        }

        # These matches have no other information
        if ( $ids[0] =~ /^(pfam|COG|cd)/ ) {
            $feat->remove_tag('product');
            $feat->add_tag_value( 'product', 'hypothetical protein' );
        }

        # Change cluster, score, and group tags to a note tag if there's data
        my @clusters = $feat->remove_tag('cluster')
          if $feat->has_tag('cluster');
        my @groups = $feat->remove_tag('group') if $feat->has_tag('group');
        my @scores = $feat->remove_tag('score') if $feat->has_tag('score');

        if ( $clusters[0] && $groups[0] && $scores[0] ) {
            my $note =
              "similar to CDD " . $clusters[0] . ", group " . $groups[0];
            $feat->add_tag_value( 'note', $note );

            #$note = "RPS BLAST score=" . $scores[0];
            #$feat->add_tag_value('note',$note);
        }
    }

    $feat;
}

sub remove_semicolon {
    my ( $product, $feat ) = @_;

    # Remove clauses with semicolons from 'product', add to a 'note'
    if ( $product =~ /^([^;]+);\s*(.+)/ ) {
        $product = $1;
        $feat->add_tag_value( 'note', $2 );
    }

    ( $product, $feat );
}

sub remove_species {
    my ( $product, $feat ) = @_;

    # If product looks something like:
    # "GpL [Enterobacteria phage P2] ..." then correct
    if ( $product =~ /(\w+)\s+\[(Enterobacteria[^]]+)\]/ ) {
        $product = "$2 $1";
    }

    # No species names allowed in product/
    if ( $product =~ /^Similar to\s+(.+?)\s+of\s+S+\s+\S+/ ) {
        my $note = $product;
        $product = $1;
        $product =~ s/^(putative|probable)\s+//;
        $feat->add_tag_value( 'note', $note );
    }

    ( $product, $feat );
}

sub fix_cog {
    my ( $product, $feat ) = @_;

    # If product name comes from match to COG, not UniRef
    if ( $product =~ /^(COG\d+):\s+(.+)/ ) {
        $product = $2;
        my $note = "similar to $1";
        $feat->add_tag_value( 'note', $note ) if ( !$1 eq 'family' );
    }

    ( $product, $feat );
}

sub fix_uniref {
    my ( $product, $feat ) = @_;

# If the product comes from UniRef, something like:
# "UPF0076 protein yjgF n=146 Tax=Bacteria RepID=YJGF_ECOL6" or
# "Putative Orf27; P2 LysB homolog; control of lysis [Ente. n=2 Tax=Yersinia RepID=Q66BL7_YERPS"
    if ( $product =~ /(.+?)\s+n=\d+\s+Tax=(.+)/ ) {

        my ( $newProduct, $species ) = ( $1, $2 );
        $species =~ s/RepID/UniRef RepID/;

        my $note = "similar to $newProduct of $species";
        $note =~ s/\sn=\d+\s/ /;

       # Remove COG or FOG ids found in the UniRef headers, e.g.:
       # FOG: TPR repeat n=1 Tax=Vibrio vulnificus RepID=Q8DF47_VIBVU
       # COG0784: FOG: CheY-like receiver n=1 Tax=Bacillus anthracis RepID=UPI00
        $newProduct =~ s/^(FOG:\s+|COG\d+:\s+FOG:\s+|COG\d+:\s*)//;

        $feat->add_tag_value( 'note', $note );

        # Add back the product but do not use locus tags from other
        # genomes, e.g. product="hypothetical protein YPO0973
        if ( $newProduct =~ /^hypothetical protein\s+.*/i ) {
            $product = 'hypothetical protein';
        }
        else {
            $product = $newProduct;
        }

        $feat->remove_tag('score') if $feat->has_tag('score');
    }

    ( $product, $feat );
}

sub correct_spelling {
    my $product = shift;

my @strs = (
'Lipoprotein_5', 'Lipoprotein 5',
'Catabolite protein activator', 'catabolite protein activator',
'Accessory protein regulator ([A-Z])', 'accessory protein regulator $1',
'Penicillin V acylase\. Cysteine peptidase\.', 'penicillin V acylase; cysteine peptidase;',
'FtsH-2 peptidase. Metallo peptidase. MEROPS family M41,', 'FtsH-2 peptidase; Metallo peptidase; MEROPS family M41;',
'Camelysin\.', 'Camelysin;',
'Metallo peptidase\.', 'Metallo peptidase;',
'D-Ala-D-Ala carboxypeptidase A\.', 'D-Ala-D-Ala carboxypeptidase A;',
'Serine peptidase\.', 'Serine peptidase;',
'carboxypeptidase DacF\.', 'carboxypeptidase DacF;',
'acetoincleaving', 'acetoin cleaving',
'dyhydrogenase', 'dehydrogenase',
'carrierprotein', 'carrier protein',
'proteinl', 'protein',
'POLY\(A\) POLYMERASE \' TRNA NUCLEOTIDYLTRANSFERASE', 'tRNA nucleotidyltransferase',
'MOSQUITOCIDAL TOXIN PROTEIN', 'mosquitocidal toxin protein',
'MONO-ADP-RIBOSYLTRANSFERASE C3', 'mono-ADP-ribosyltransferase C3',
'hyothetical', 'hypothetical',
'biosynthsis', 'biosynthesis',
'AMIDOHYDROLASE', 'amidohydrolase',
'ACETOIN TRANSPORT PERMEASE PROTEIN', 'acetoin transport permease protein',
'(trasporter|transoprted|tranporter)', 'transporter',
'3-OXOADIPATE ENOL-LACTONASE', '3-oxoadipate enol-lactonase',
'4\'\'\'\'-phosphopantetheinyl', 'phosphopantetheinyl',
' \(And other\) ', ' ',
'sensor\(S\)', 'sensors',
'TRANSPOSON', 'Transposon',
'INVERTASE', 'Invertase',
'\bsopre\b', 'spore',
'proteintic', 'protein',
'NAD\(P\)H', 'NADPH',
'spaning', 'spanning',
'proetin', 'protein',
'Hypotethical', 'hypothetical',
'reguatory', 'regulatory',
'fibre', 'fiber',
'Uncharacterised', 'uncharacterized',
'putaive', 'putative',
'haemin', 'hemin',
'Haemolytic', 'Hemolytic',
'unknow ', 'unknown ',
'haemagglutinin', 'hemagglutinin',
'PERIPLASMIC PROTEIN', 'periplasmic protein',
'MITOCHONDRIAL TRANSPORTER ATM1 \(Atm1\)', 'Mitochondrial transporter Atm1',
'Protein C\. Serine peptidase\.', 'Protein C serine peptidase',
'unkown', 'unknown',
'addtional', 'additional',
'protein protein', 'protein',
'similar to Similar', 'similar',
'^Acyl\stransferasesregion$', 'Acyl transferase',
'permease-associated\sregion$', 'permease domain protein',
'Fragment', 'fragment',
'Putative', 'putative',
'characteris', 'characteriz',
'[Uu]ncharacterized', 'putative',
'[Ss]ulphate', 'sulfate',
'- \(pentapeptide', '-(pentapeptide',
'genes activator', 'gene activator'
);

    while ( my ($search,$replace) = splice(@strs,0,2) ) {
        $product =~ s/$search/$replace/;
    }

    $product;
}

sub remove_loci {
    my $product = shift;

    my @strs = (
'(Uncharacterized adenine-specific methylase)\s+[\d{1,}a-z{1,}\/-_]+',
'(hydrolase)\s+[\d{1,}a-z{1,}\/-_]+',
'(metalloprotease)\s+[\d{1,}a-z{1,}\/-_]+',
'(Phosphodiesterase),?\s+[\d{1,}a-z{1,}\/-_]+',
'(reductase)\s+[\d{1,}a-z{1,}\/-_]+',
'(transport protein)\s+[\d{1,}a-z{1,}\/-_]+',
'(phosphotransferase)\s+[\d{1,}a-z{1,}\/-_]+',
'(Uncharacterised conserved protein)\s+[\d{1,}a-z{1,}\/-_]+',
'(Uncharacterized MFS-type transporter)\s+[\d{1,}a-z{1,}\/-_]+',
'(Uncharacterized transporter)\s+[\d{1,}a-z{1,}\/-_]+',
'(Pirin-like protein)\s+[\d{1,}a-z{1,}\/-_]+',
'(Uncharacterized protease)\s+[\d{1,}a-z{1,}\/-_]+',
'(HAD-superfamily hydrolase, subfamily IB,)\s+[\d{1,}a-z{1,}\/-_]+',
'(Uncharacterized mscS family protein)\s+[\d{1,}a-z{1,}\/-_]+',
'(Shikimate 5-dehydrogenase-like protein)\s+[\d{1,}a-z{1,}\/-_]+',
'(Peptidase T-like protein)\s+[\d{1,}a-z{1,}\/-_]+',
'(Uncharacterized ABC transporter ATP-binding protein)\s+[\d{1,}a-z{1,}\/-_]+',
'(Uncharacterized acyl-CoA thioester hydrolase)\s+[\d{1,}a-z{1,}\/-_]+',
'UPF\d+\s+(\w+-binding protein)',
'^([^()]+ \([^)]+\)) \([^)]+\)'
    );

    # No loci in product names, e.g 'hydrolase LC_123'
    for my $str (@strs) {
        return $1 if ( $product =~ /$str/i );
    }

    $product;
}

sub add_trailing {
    my $product = shift;

    my @strs = (
'^\s*\w+\s+repeat\s*$',                       '-containing protein',
'transcriptional\sregulator,\sPadR-like\s*',  ' protein',	
'DinB\sfamily\ssuperfamily\s*',               ' protein',
'YfhF-\s*',                                   '-like',
'WbqC-\s*',                                   '-like',
'Zinc\sfinger,\sCHC2-\s*$',                   '-like',
'Tetratricopeptide\sTPR_4\s*$',               ' containing protein',
'sensor\shistidine\skinase\sdomain\s*$',      ' containing protein',
'SWIM\szinc\sfinger\s*$',                     ' containing protein',
'Tetratricopeptides\(TPR\)\srepeat\s*$',      ' containing protein',
'transporter\srelated$',                      ' protein',
'^(\w+)\sbinding\sdomain$',                   ' protein',
'(\S+-like)$',                                ' protein',
'SpoVR like',                                 ' protein',
'Glutathione S-transferase domain',           ' protein'
);

    while ( my ($str,$add) = splice(@strs,0,2) ) {
        $product .= $add if ( $product =~ /$str/i );
    }

    $product;
}



sub remove_similar {
	my $product = shift;

my @strs = (
'^Similar to (acetyltransferase|exochitinase|restin isoform b|permease protein of ABC transport system|protein gp49 from prophage N15)',
'^Similar to (fimbrial subunit type 1|bacteriophage integrase|vgrG protein|base plate protein gp25 of Bacteriophage|bacteriophage tail fiber assembly protein|protein V)',
'^Related to (galactoside O-acetyltransferase)/',
'^Strongly similar to (S-adenosylmethionine synthetase)',
'^Exhibited considerable similarity to a small (polypeptide \(RepA\))',
'^Similarities with (transcription regulator LysR family|alpha replication protein of prophage|Photorhabdus cytotoxin|tail fiber protein)'
);

	# Remove phrases like 'similar to...'
    for my $str (@strs) {
        return $1 if $product =~ /$str/i;
    }

	$product;
}

sub make_singular {
	my $product = shift;

	my @strs = (
'GTPases', 'GTPase',
'variants', 'variant',
'synthetases', 'synthetase',
'glucosidases', 'glucosidase',
'thioredoxins', 'thioredoxin',
'asparaginases', 'asparaginase',
'acetylases', 'acetylase',
'enzymes', 'enzyme',
'Flavodoxins', 'Flavodoxin',
'toxins', 'toxin',
'Permeases', 'Permease',
'components', 'component',
'proteins', 'protein',
'systems', 'system',
'regulators', 'regulator',
'phosphatases', 'phosphatase',
'determinants', 'determinant',
'recombinases', 'recombinase',
'transferases', 'transferase',
'ATPases', 'ATPase',
'exporters', 'exporter',
'hydrolases', 'hydrolase',
'reductases', 'reductase',
'cytochromes', 'cytochrome',
'proteases', 'protease',
'kinases', 'kinase',
'transporters', 'transporter',
'oxidases', 'oxidase',
'helicases', 'helicase',
'synthases', 'synthase',
'peptidases', 'peptidase',
'Dehydrogenases', 'Dehydrogenase',
'lyase.+?and.+?lyases', 'lyase',
'protein protein', 'protein',
'solvents', 'solvent'
);

    while ( my ($search,$replace) = splice(@strs,0,2) ) {
        $product =~ s/$search/$replace/ig;
    }

	$product;
}

sub edit_note {
    my $feat = shift;

    if ( $feat->has_tag('note') ) {
        my @notes = $feat->remove_tag('note');
        for my $note ( @notes ){
            next if ( $note =~ /hypothetical protein/ );
            $note = remove_banned($note);
            $feat->add_tag_value('note',$note);
        }
    }     

    $feat;
}

sub remove_banned {
	my $product = shift;

	my @strs = (
'\s+related\s+',            
'\s+homologs?\s*',                           
'\s+\(partial\s\)',                          
',?\s*putative',                          
'^(Probable|Possible|Predicted|Putative)\s+',
'\s+\(Fragment\)\s?',
'\bgene\b/\bprotein\b',                      
'^PREDICTED:\s*',                    
'^(B.thurinienis|Salmonella)\s+',  
'^Similar to\s+',           
'^Truncated\s+',
'hypothetical protein'                      
	);

    for my $str ( @strs ) {
        $product =~ s/$str/ /ig;
    }
    $product = trim($product);

	$product;
}

sub remove_trailing {
    my $product = shift;

    # Remove trailing non-text
    $product = trim($product);
    $product =~ s/(_|,|\?)$//;

    # Remove meaningless trailing comments
    my @strs = (
', PFL_4704',
' of prophage CP-933K',
',\s+PFL_\d+',
'\s+and\sother\spencillin-binding\sproteins?',
'\s+\(Insoluble\sfraction\)',
'\s+\(amino terminus\)',
'\s+(SA2311|A118|AAA_5|MJ\d+|YLL\d+|SAB\d+[a-z]+|alr\d+)',
'\/Dioxygenas',
's+domain\s1\s\(GGDEF\)',
'\s+C\sX\sregion',
'\s+\(subunit\)',
'\s+firmicutes',
'\(Dihydro>',
'\s+in\S+\s+\S+region',
'\s+catalytic region',
'\s+and\s+(dehydrogenase|aminotransferase)',
's+\(Glutamate-aspartate\scarrierprotein\)',
'\s+\(Replication\sprotein\sori\d+\)',
'\s+\(Replication\sprotein\)',
'in\sMarinococcus\shalophilus',
'\s+clpC\/mecB',
'\s+CA_[A-Z]+\d+',
',?\s+fhuD',
'\s+BCE_\d+',
'\s+BLi\d+\/BL\d+',
'\s+(BT|LMOf)\d+_\d+',
'family\sprotein',
'HD\ssub\sdomain',
'\s*(SA|BH|VC|NMB|HI|SH)\d+',
'\s+pXO\d+-\d+\/BXB\d+\/GBAA_pXO\d+_\d+',
'\s+BA_\d+\/GBAA\d+\/BAS\d+',
'\s*\(Ans\soperon\srepressor\s*protein\)',
'\s*\(Divided\swith\sOB2865\sand\sOB2866\)',
'\s*\(N-terminal\sdomain\sto\sN-Acetylmuramoyl-L-alanine\samidase\)',
'\s+family\sfamily',
'\s+related',
',?\sfamily',
'\s*of\sV.\sanguillarum\s\(YclQ\sprotein\)',
'\s*\(Putative\sendopeptidase\sinhibitor\)',
'\s*\(Putative\sarconitate\shydratase\)',
'\s*\(Two component\ssystem\sresponse\sregulatory\sprotein\)',
'\s*\(Hypothetical\sprotein\)',
'\s*\(Outer\smembrane\susher\sprotein\)',
'\s*\(eIF-2Bgamma.eIF-2Bepsilon\)',
'\s*\(some\scontain\sLysM.invasin\sdomains\)',
',\scatalytic\sdomain:D-isomer\sspecific\s2-hydroxyacid\sdehydrogenase,\sNAD\sbinding\sdomain',
',?\struncat(ion|ed)',
',?\sinterruption',
',?\sYHCS\sB.subtilis\sortholog',
',?\s\([A-Z\d-]+\)',
',\ssimilar\sto\sSW:\w+',
',?\shexapeptide\srepeat',
',?\s(C|N)-termin(us|al)',
',?\s(N|C)-terminal\s(domain|region)',
',?\s(SHL|HD)\sdomains?',
'HD\ssub\sdomain',
',?\s+N-region',
'and\sinactivated\sderivatives',
'\s+and\s+orf\d+',
'\s+and',
',?\s(mitochondrial|chloroplastic)',
'\s*BCE?_\d+',
'\(Dihydrolipoamide\sacetyltransferase\scomponent\sof\sacetoin\scleavingsystem\)',
'\s*\(Acetoin\sdehydrogenase\sE2\scomponent\)',
'\s*,\ssmall\schain\sVC\d+',
'\s*,\sGNAT\sfamily\sfamily\sprotein',
'\s*BC_\d+',
'\s*,\speriplsmic\sligand-binding\sprotein',
'\s*BA_\d+\/GBAA\d+\/BAS\d+',
'\s*(\[NAD\(P\)+\]|\[Mn\]|\[NAD\(P\)H\]|\[NADPH?\]|\[ATP\]|\[Cu-Zn\]|\[ATP hydrolyzing\]|\[isomerizing\]|\[carboxylating\]|\[glutamine-hydrolyzing\]|\[decarboxylating\]|\[a?symmetrical\])',
'\/isomerase:Polysaccharide biosynthesis protein CapD:dTDP-4-dehydrorhamnose reductase:Nucleotide sugar epimerase',
'\s+\(Diaminohydroxyphosphoribosylaminopyrimidine deaminase \(Riboflavin-specific deaminase\) and 5-amino-6-\(5-phosphoribosylamino\)uracil reductase\)',
'\s+\(Probable\), IS891\/IS1136\/IS1341:Transposase, IS605 OrfB',
'\s+and inactivated derivatives-like protein',
'\[cytochrome\]\s*',
', fused inner membrane subunits',
', auxiliary component',
', periplasmic component',
', transcription of rRNA and tRNA operons, and DNA replication',
' HI_\d+',
'\/FOG: TPR repeat protein',
'\/RND superfamily resistance-nodulation-cell division antiporter',
'\/ DNA internalization-related competence protein ComEC\/Rec2 \/ predicted membrane metal-binding protein',
', transcription of rRNA and tRNA operons, and DNA replication',
': membrane component/ATP-binding',
' MS\d+',
', HI\d+',        
' HD_\d+',
', regulator of competence-specific genes',
', fused lipid transporter subunits of ABC superfamily.+'
);

    for my $str (@strs) {
        $product =~ s/$str$//i;
    }

    $product;
}

sub is_hypothetical {
	my $product = shift;

	return 1 if ! $product;

	# Product is 'hypothetical' if just an id with a vague name, 
	# e.g. 'Lin0391 protein' ^Lin\d+ protein \(Lin\d+ protein\)$/ 

	my @strs = (
'^RHS$',
'^protein HIB_\d+$',
'^Ybl\d+$',
'^(similar|similarities)\s+(to|with)\s+(unknown|putative|probable|C-terminal|N-terminal)',
'polypeptide \([a-z]+\)\sencoded\sby\splasmid',
'^Possible\speptide\santibiotic',
'^PugilistDominant',
'mutants\sblock\ssporulation\safter\sengulfment',
'^Homo\ssapiens',
'^Similar\sto\sORF13\sof\senterococcus\sfaecalis',
'^ORF13\sof\senterococcus\sfaecalis\sTRANSPOSON\sTN916',
'^CheY-homologous\sreceiver\sdomain',
'^Plasmid\spPOD2000\sRep,\sRapAB,\sRapA,\sParA,\sParB,\sand\sParC',
'^Plasmid\spRiA4b\sORF-3',
'DEHA0C09658g\sDebaryomyces\shansenii',
'Bacillus\scereus\sgroup-specific\sprotein',
'^UPF\d+.+?protein.+?\S+$',
'^BC\d+\w+\sprotein',
'N\sterminal\sregion\sof\sphage-related\scontractile\stail\ssheath\sprotein',
'chromosome\s+\d+open\s+reading\s+frame\s+\d+',
'complete\sgenome',
'Genome\ssequencing\sdata,\scontig\sC\d+',
'chromosome\s\d+\sopen\sreading\sframe\s\d+',
'complete\snucleotide\ssequence',
'DNA,\scomplete\ssequence',
'^Genomic\sDNA',
'whole\sgenome\sshotgun\ssequence',
'Gene,\scomplete\scds',
'^Orf\s+\d+[A-Z]+$\s',
'^Nucleoside\srecognition',
'^Required\sfor\splasmid\sstability$',
'^Possible\sBacterial\sIg-like\sdomain',
'^Alpha-peptide$',
'LPXTG-motif\scell\swall\sanchor\sdomain$',
'^Amino\sacid\stranporter$',
'^Ankyrin\srepeats\scontaining\sprotein$',
'^Biotin\lipoyl\sattachment$',
'^Antigen$',
'^Restriction\smodification\ssystem\sDNA\sspecificity\sdomain',
'^ABC\s\(ATP-binding\scassette\)\stransporter\snucleotide-binding\sdomain$',
'^(SAF|NACHT|Resolvase|FRG|C1q)\s+domain$',
'^Contains\scell\sadhesion\sdomain$',
'^Gene,\sIS3-like\selement$',
'^Micrococcal\snuclease-like\sprotein$',
'^modification\smethylase\sOB\d+$',
'^Micrococcal\snuclease$',
'^leucine-rich\srepeat-containing\sprotein\sDDB\d+$',
'^Possible\ssensor\shistidine\skinase\sdomain',
'^PIN\s.PilT\sN\sterminus.\sdomain$',
'^Divergent\sAAA\sregion$',
'^Conserved\srepeat\sdomain\sprotein$',
'^Collagen\striple\shelix\srepeat$\s',
'^Parallel\sbeta-helix\srepeat$\s\s',
'^PBS\slyase\sHEAT-like\srepeat\s\s',
'^Sel1-like\srepeat$',
'^Parallel\sbeta-helix\srepeat$\s\s',
'^Helix-turn-helix\smotif$',
'^Transferase\shexapeptide\srepeat$',
'^Helix-turn-helix,\stype\s\d+$',
'^(Thioredoxin|HTH|Potential\sSulfotransferase)\sdomain$',
'^S-layer\sdomain$',
'^(Thioredoxin|HTH)\sdomain\sfamily$',
'^Amino\sacid\sadenylation\sdomain$\s',
'^CopG-like\sDNA-binding$',
'^Helix-turn-helix$',
'^Helix-turn-helix\sHxlR\stype$',
'^Alpha\beta\shydrolase\sfold$\s\s',
'^Hydrolase,\salpha.beta\sfold\s\s',
'^Hydrolase,\salpha.beta\shydrolase\sfold',
'^unknown\s\w+-like\sprotein$',
'^C-terminal\shalf\sof\sCry\sprotein$\s',
'^protein\sof\sunknown\sfunction\s',
'^unknown\sprotein\s[A-Z]\d+\s\s',
'^Similarities\sto\sphage\sprotein\s',
'^hypothetical\smembrane\sspanning\sprotein$\s',
'^Conserved\shypothetical\sprotein$\s',
'^Uncharacterized\sconserved\smembrane\sprotein\s',
'^Similar\sto\s(hypothetical|bacteriophage)\sprotein\s',
'^Similar\sto\sshort\sunknown\sprotein\s',
'^Truncated\sphage-like\s',
'^Hypothetical\sphage\sprotein$\s',
'^Phage-related\sprotein\s',
'^Hypothetical\stranscriptional\sregulator$\s',
'^Mitochondrial\stransporter\s',
'^Uncharacterized\slow-complexity\s',
'^Conserved\shypothetical\smembrane\sprotein$\s',
'^integral\smembrane\sprotein\sTIGR\d+$\s',
'^hypothetical\stwo\sdomain\sprotein$\s',
'^Conserved\shypothetical\sintegral\smembrane\sprotein\s',
'^[a-z]+\s\(Conserved\sprotein\s[a-z]+\)$\s',
'^[a-z]+\s\(Conserved\smembrane\sprotein\s[a-z]+\)$\s',
'^Membrane\sprotein,\sputative$\s',
'^Conserved\sprotein\s*$\s',
'Predicted\stranscriptional\sregulator\swith\san\saddtional\sconserved\sdomain\s',
'transcriptional\sregulator\swith\san\sadditional\sconserved\sdomain\s',
'^Conserved\s(predicted|phage|domain|membrane|exported|hypothetical)\sprotein\s',
'^membrane\sspann?ing\sprotein$\s',
'^inner\smembrane\sprotein$\s',
'^Predicted\ssmall\ssecreted\sprotein$\s',
'^Uncharacterized\sconserved\ssmall\sprotein\s',
'^Uncharacterized\sconserved\sprotein\s(UCP|CAC)\d+\s',
'^Outer\smembrane\sprotein\s[A-Z]+\d+$\s',
'^Predicted\sprotein\s*$\s',
'^Uncharacterized\sconserved$\s',
'putative\suncharacterized\sprotein\s',
'^Hypothetical\s(conserved|exported)\sprotein$\s',
'^Hypothetical\sprotein,\spartial$\s',
'^hypothetical\sprotein\smembrane\sprotein$\s',
'^Uncharacterized\sprotein\sconserved\sin\sbacteria$\s',
'^Uncharacterized\sconserved$\s',
'^lipoprotein$\s',
'^antigen$\s',
'^Bifunctional\senzyme,\scontains\s',
'^Group-specific\sprotein$\s',
'^Extended\sORF\sof\s',
'^Uncharacterized\s(membrane|conserved)\sprotein$\s',
'^Exported\smembrane\sprotein$\s',
'^Uncharacterized\sinner\smembrane\sprotein$\s',
'^Uncharacterized\sconserved\ssmall\sprotein-like\sprotein\s',
'^Predicted\stranscriptional\sregulator\swith\san\saddtional\sconserved\sdomain\s',
'^Uncharacterized\sHTH-type\stranscriptional\sregulator\s',
'^Uncharacterized\sprotein\s',
'^uncharacterized\sdomain\s\d+$\s',
'^Uncharacterized\s[.\d]+\skDa\sprotein\s',
'^Possible\s(CDF|PET)\sfamily\s',
'^(Possible|Similarity|Uncharacterized|protein|toxin|predicted|hypothetical|exported)$\s',
'^\w+-related$\s',
'^Enzyme\s?$\s',
'^Transposase\sIS\d+.$\s',
'^Bacillus\scereus\sspecific\sprotein,\suncharacterized$\s',
'^[A-Z][a-z]{2}[A-Z]$\s',
'^LP\d+G.\d+$\s',
'^R.EcoHK31I\sprotein$\s',
'^lmo\d+$\s',
'^Orf\d+\s*$\s',
'^\d+orf\d+$\s',
'^Orf\s+\d+[A-Z]+\s*$\s',
'^Orf\s+[A-Z]+\d+\s*$\s',
'^(IS66\sOrf2|Orf2)\slike$\s',
'^Orf\s+[A-Z]+\d+\s+putative\s',
'^IS66\sOrf2\slike$\s',
'^Orf\d+-like\sprotein\s*$\s',
'^Orf\d+-[A-Z]+\d+\s*$\s',
'^Orf\S+\s+\(\S+\s+protein\)\s*$\s',
'^Ig\s+hypothetical\s+\d+$\s',
'^[a-z]{2}$\s',
'^(HI|EF|Orf|Mob|MW|Blr|Cro|p|Blr|Orf-)\d+;?$\s',
'^Possible\s\(\S+\)\sorf\d+\s',
'^EF\d+\s+\(.+?\)$\s',
'^Orf[A-Z]\s*$\s',
'^GG\d+,?\s+isoform\s+[a-z]+$\s',
'^\d+orf\d+\s*$\s',
'^Cj\d+-\d+$\s',
'^PXO\d+-\d+\s*$',
'^UPI\S+\s+cluster$\s',
'^Vng\d+[a-z]$\s',
'^Zwa5b$\s',
'^PclL$\s',
'^\d+$\s',
'^(Rep1|Doc|Vip2Ac|Zwa5A|Ugd|Sip1A|Vip1A\(BR\)|Tnp166|Blr5358|LAAC|MW2125)$',
'^Gene\s\d+\sprotein$\s',
'^Gp\d+\.\d+$\s',
'^(19|Gne|Aec\d{1,2})$\s',
'^protein\s+[\d{1,}a-z{1,}\/_-]+$\s',
'^[\d{1,}a-z{1,}\/_-]+\s+protein$\s',
'^lipoprotein\s+[\d{1,}a-z{1,}\/-_]+$\s',
'^Maf-like\sprotein\s+[\d{1,}a-z{1,}\/-_]+$\s',
'^[\d{1,}a-z{1,}\/-_]+\s+protein\s+[\d{1,}a-z{1,}\/-_]+$\s',
'Eag\d+',
'Conserved domain protein',
'Hypothetical conserved protein',
'Uncharacterized conserved protein',
'Uncharacterized membrane protein',
'^protein$',
'^gene$',
'^genes$',
'^conserved protein$',
'UPI[\dA-Z]+\s+cluster',
'KLLA0D19929p',
'protein PARA_\d+',
'HI\d+-like protein',
'protein conserved in bacteria',
'protein PM\d+',
'membrane protein PM\d+',
'protein HI_\d+',
'Hybrid protein containing carboxymuconolactone decarboxylase domain'
);

	for my $str ( @strs ) {
		return 1 if ( $product =~ /$str/i );
	}

	0;
}

sub set_feat_from_adj_contig {
	my ($self,$feat) = @_;

	if ( $feat ) {
		push @{$self->{feat_from_adj_contig}},$feat;
	} else {
		$self->{feat_from_adj_contig} = undef;
	}
}

sub get_feat_from_adj_contig {
	my $self = shift;

	if ( defined $self->{feat_from_adj_contig} ) {
		return @{$self->{feat_from_adj_contig}};
	}

	return ();
}

sub get_start_plus {
	my ($self,$geneend) = @_;

	my $mod = $geneend % 3;

	$mod + 1;
}

sub get_start_minus {
	my ($self,$fend,$cend) = @_;
	my $diff = $fend - $cend;

	return 1 if ( ($diff % 3) == 0 ); # Not tested
	return 2 if ( ($diff % 3) == 2 ); # Not tested
	return 3 if ( ($diff % 3) == 1 ); # Confirmed
}

sub newFeatures {
	my ($self, @feat) = @_;

	push @{$self->{newFeatures}}, @feat if @feat;
	return @{$self->{newFeatures}};
}

# Example command, 4/2012:
# tbl2asn -t template -p path_to_files -X C -M n -Z discrep \
# -j "[organism=Clostridium difficile ABDC] [strain=ABDC] [host=Homo sapiens] [country=Canada: Montreal] \
# [collection_date=2008] [isolation-source=stool sample] [note=isolated from an outbreak in Montreal] [gcode=11]"
sub run_tbl2asn {
	my ($self,$comment,$run) = @_;

    my ( $tmplt, $outdir, $tbl2asn, $id, $gcode, $strain, $organism, $host,
        $country, $collection_date, $isolation_source, $submission_note )
      = (
        $self->template,         $self->outdir,
        $self->executable,       $self->id,
        $self->gcode,            $self->strain,
        $self->organism,         $self->host,
        $self->country,          $self->collection_date,
        $self->isolation_source, $self->submission_note
      );

	if ( $run > 1) {
		for my $suffix ( qw( val sqn gbf ) ) {
			system "mv $outdir/$id.$suffix $outdir/$id.$suffix.orig" if ( -e "$outdir/$id.$suffix" );
		}
		system "mv discrp discrp.orig" if ( -e "discrp" );
	}

  # Fix country and city
  $country =~ s/:\s*/: /;

  my $jstring = "[organism=$organism $strain] [strain=$strain] [gcode=$gcode]";
  $jstring   .= " [host=$host]" if $host;
  $jstring   .= " [country=$country]" if $country;
  $jstring   .= " [collection_date=$collection_date]" if $collection_date;
  $jstring   .= " [isolation-source=$isolation_source]" if $isolation_source;
  $jstring   .= " [note=$submission_note]" if $submission_note;

	my $cmd = "$tbl2asn -t $tmplt.sbt -p $outdir -M n -Z discrp -y \"$comment\" " .
                  "-X C -V b -w $outdir/$id.cmt -j \"$jstring\"";
	print "tbl2asn command is \'$cmd\'\n" if $self->debug;
	`$cmd`;
	return 1;
	
	die "Problem running tbl2asn: $cmd";
}


sub fix_discrp {
	my $self = shift;

	my $tbl = $self->read_tbl;
	print "Fixing discrepancies\n" if $self->debug;

	my @geneoverlaps = $self->get_gene_overlaps;
	$tbl = $self->delete_from_tbl($tbl,@geneoverlaps);

	my @rnaoverlaps = $self->get_rna_overlaps;
	$tbl = $self->delete_from_tbl($tbl,@rnaoverlaps);

	my @duprnas = $self->get_dup_rnas;
	$tbl = $self->delete_from_tbl($tbl,@duprnas);

	$self->write_tbl($tbl);

	1;
}

# DiscRep:OVERLAPPING_GENES::2 genes overlap another gene on the same strand.
# Gene    ykris_r40       lcl|contig00434:49-1579 ykris_r40
# Gene    ykris_r30       lcl|contig00434:1101->1648      ykris_r30
# or:
# DiscRep_ALL:FIND_OVERLAPPED_GENES::6 genes completely overlapped by other genes
# WGQ:Gene        WGQ_50  lcl|ctg7180000000008:7373-7807  WGQ_50
# WGQ:Gene        WGQ_240 lcl|ctg7180000000008:33916-34350        WGQ_240
sub get_dup_rnas {
	my $self = shift;
	my @todelete = ();
	my $id = $self->id;
	my $gene1;

	my @overlapgenes   = $self->get_from_discrp('OVERLAPPING_GENES');
	# push @overlapgenes, $self->get_from_discrp('FIND_OVERLAPPED_GENES');
	print "Overlapping genes: @overlapgenes\n" if $self->debug;

	$gene1 = shift @overlapgenes;

	while ( my $gene2 = shift @overlapgenes ) {

		if ( $gene1 =~ /_r\d+/ && $gene2 =~ /_r\d+/ ) {

			my ($g1start,$g1end) = $gene1 =~ /(?:contig|ctg)[\d.]+:c?(\d+)[<>]?-[<>]?(\d+)/;
			my ($g2start,$g2end) = $gene2 =~ /(?:contig|ctg)[\d.]+:c?(\d+)[<>]?-[<>]?(\d+)/;

			my ($g1ctg) = $gene1 =~ /((?:contig|ctg)[.\d]+)/;
			my ($g2ctg) = $gene2 =~ /((?:contig|ctg)[\d.]+)/;

			# Remove identical rRNAs and rRNAs inside rRNAs
			if ( $g1start ==  $g2start && $g1end == $g2end && $g1ctg eq $g2ctg ) {
 				push @todelete, $gene2;
 				print "Will remove duplicate gene2: $gene2" if $self->debug;
			} elsif ( $g1start >=  $g2start && $g1end <= $g2end && $g1ctg eq $g2ctg ) {
				push @todelete, $gene1;
				print "Will remove smaller rRNA inside larger: $gene1" if $self->debug;
			} elsif ( $g2start >=  $g1start && $g2end <= $g1end && $g1ctg eq $g2ctg ) {
		 		push @todelete, $gene2;
				print "Will remove smaller rRNA inside larger: $gene2" if $self->debug;
			}

			my $g1len = abs($g1start - $g1end);
			my $g2len = abs($g2start - $g2end);

			# Remove the smaller of 2 overlapping rRNAs
			if ( $g1len >= $g2len ) {
		 		push @todelete, $gene2;
				print "Will remove smaller rRNA overlapping larger: $gene2" if $self->debug;
			} else {
				push @todelete, $gene1;
				print "Will remove smaller rRNA overlapping larger: $gene1" if $self->debug;
			}

		}

		$gene1 = $gene2;
	}
	@todelete = unique(@todelete);
	@todelete;
}

# "DiscRep:OVERLAPPING_GENES::277 genes overlap another gene on the same strand."
# "DiscRep:OVERLAPPING_CDS::12 coding regions overlap another coding region with a 
# similar or identical name."
# "DiscRep:OVERLAPPING_CDS::8 coding regions overlap another coding region with a similar 
# or identical name that does not contain 'ABC' and do not have the appropriate note text"
# "DiscRep:OVERLAPPING_CDS::4 coding regions overlap another coding region with a similar 
# "or identical name that contains 'ABC'"
# Solution: where 1 gene completely overlaps another, delete the smaller gene.
# Note that the DiscRep:OVERLAPPING_GENES class contains all the smaller classes
sub get_gene_overlaps {
	my $self = shift;
	my @todelete = ();
	my $id = $self->id;
	my $gene1;

	my @overlapgenes = $self->get_from_discrp('OVERLAPPING_GENES');
	my @overlapcds = $self->get_from_discrp('OVERLAPPING_CDS');
	push @overlapgenes,@overlapcds;

	$gene1 = shift @overlapgenes;

# 192:CDS    glycosyltransferase lcl|ctg7180000000003:c863978-863178 192_7700
# 192:CDS    glycosyltransferase lcl|ctg7180000000003:c864862-863978 192_7710

	while ( my $gene2 = shift @overlapgenes ) {

		# Note that this loop does not collect rRNAs or tRNAs

		my ($g1start,$g1end) = $gene1 =~ /(?:contig|ctg)[\d.]+:c?(\d+)[<>]?-[<>]?(\d+)/;
		my ($g2start,$g2end) = $gene2 =~ /(?:contig|ctg)[\d.]+:c?(\d+)[<>]?-[<>]?(\d+)/;

		my $g1len = abs($g1start - $g1end);
		my $g2len = abs($g2start - $g2end);

		my ($g1ctg) = $gene1 =~ /((?:contig|ctg)[.\d]+)/;
		my ($g2ctg) = $gene2 =~ /((?:contig|ctg)[\d.]+)/;

		# Gene	yberc_40220	lcl|contig01136:c856-281	yberc_40220
		# Gene	yberc_40230	lcl|contig01136:c1053-856	yberc_40230
 		if ( $g1start >=  $g2end && $g1start <= $g2start && $g1ctg eq $g2ctg ) {
			# Delete the small one
 			if ( $g1len >= $g2len ) {
 				push @todelete,$gene2 if ( $gene2 !~ /${id}_(r|t)/ );
 				print "1a: Will remove gene2: $gene2" if $self->debug;
 			} else {
 				push @todelete,$gene1 if ( $gene1 !~ /${id}_(r|t)/ );
 				print "1b: Will remove gene1: $gene1" if $self->debug;
 			}
 		}

		# Gene	yberc_39260	lcl|contig00890:6892-7563	yberc_39260
		# Gene	yberc_39270	lcl|contig00890:7563-7889	yberc_39270
 		if ( $g1end >=  $g2start && $g1end <= $g2end && $g1ctg eq $g2ctg ) {
			# Delete the small one
			if ( $g1len >= $g2len ) {
				if ( $gene2 !~ /${id}_(r|t)/ ) {
					push @todelete,$gene2;
					print "2a: Will remove gene2: $gene2" if $self->debug;
				}
			} else {
				if ( $gene1 !~ /${id}_(r|t)/ ) {
					push @todelete,$gene1;
					print "2b: Will remove gene1: $gene1" if $self->debug;
				}
			}
 		}

		$gene1 = $gene2;
	}

	@todelete = unique(@todelete);
	@todelete;
}

# "2 coding regions completely contain RNAs"
# Solution: where a gene overlaps a tRNA or mRNA delete the gene.  
# DiscRep_SUB:RNA_CDS_OVERLAP::1 coding regions are completely contained in RNAs
# bcere0010:CDS hypothetical proteinlcl|contig00145:2461-2892 bcere0010_53360
# bcere0010:rRNA 23S ribosomal RNAlcl|contig00145:9-2928 bcere0010_r50
#
# DiscRep_SUB:RNA_CDS_OVERLAP::5 coding regions completely contain RNAs
# bcere0010:CDS hypothetical proteinlcl|contig00078:84-557 bcere0010_4450
# bcere0010:rRNA 5S ribosomal RNAlcl|contig00078:89-203 bcere0010_r40
#
# DiscRep_ALL:RNA_CDS_OVERLAP::7 coding regions overlap RNA features
# FATAL: DiscRep_SUB:RNA_CDS_OVERLAP::7 coding regions are completely contained in RNAs
# WGQ:CDS Cell wall-associated hydrolase  lcl|ctg7180000000008:7373-7807  WGQ_50
# WGQ:rRNA        23S ribosomal RNA       lcl|ctg7180000000008:7227-10119 WGQ_r140
# WGQ:CDS Cell wall-associated hydrolase  lcl|ctg7180000000008:33916-34350        WGQ_240
# WGQ:rRNA        23S ribosomal RNA       lcl|ctg7180000000008:33770-36662        WGQ_r100
#
# FATAL: DiscRep_SUB:RNA_CDS_OVERLAP::6 coding regions are completely contained in RNAs
# WGQ:CDS Cell wall-associated hydrolase  lcl|scf7180000000008:c220896-220462     WGQ_1940
# WGQ:rRNA        23S ribosomal RNA       lcl|scf7180000000008:218150-221042      WGQ_r90

sub get_rna_overlaps {
	my $self = shift;
	my @todelete = ();
	my $gene1;
	my $id = $self->id;

	my @overlaps = $self->get_from_discrp('RNA_CDS_OVERLAP');
	print "RNA-CDS overlaps:@overlaps\n" if $self->debug;

	$gene1 = shift @overlaps;

	while ( my $gene2 = shift @overlaps ) {

      my ($g1start,$g1end) = $gene1 =~ /(?:contig|ctg|scf)[\d.]+:c?(\d+)[<>]?-[<>]?(\d+)/;
      my ($g2start,$g2end) = $gene2 =~ /(?:contig|ctg|scf)[\d.]+:c?(\d+)[<>]?-[<>]?(\d+)/;
      my ($g1ctg) = $gene1 =~ /((?:contig|ctg|scf)[\d.]+)/;
      my ($g2ctg) = $gene2 =~ /((?:contig|ctg|scf)[\d.]+)/;

		if ( $g1start >=  $g2end && $g1end <= $g2start && $g1ctg eq $g2ctg 
			   && $gene1 =~ /^$id:(t|r)RNA/ && $gene2 =~ /^$id:CDS/ ) {
 				push @todelete,$gene2;
 				print "RNA inside CDS, will remove CDS: $gene2" if $self->debug;
 		}
 		if ( $g2start >=  $g1end && $g2end <= $g1start && $g1ctg eq $g2ctg 
			   && $gene2 =~ /^$id:(t|r)RNA/ && $gene1 =~ /^$id:CDS/ ) {
 				push @todelete,$gene1;
 				print "RNA inside CDS, will remove CDS: $gene1" if $self->debug;
 		}
		# RNA in CDS, both +1 strand
		if ( $g2start >=  $g1start && $g2end <= $g1end && $g1ctg eq $g2ctg 
			   && $gene2 =~ /^$id:(t|r)RNA/ && $gene1 =~ /^$id:CDS/ ) {
 				push @todelete,$gene1;
 				print "RNA inside CDS, +1, will remove CDS: $gene1" if $self->debug;
 		}
 		if ( $g1start >=  $g2start && $g1end <= $g2end && $g1ctg eq $g2ctg 
			   && $gene1 =~ /^$id:(t|r)RNA/ && $gene2 =~ /^$id:CDS/ ) {
 				push @todelete,$gene2;
 				print "RNA inside CDS, +1, will remove CDS: $gene2" if $self->debug;
 		}

		# DiscRep:RNA_CDS_OVERLAP::4 coding regions overlap RNA features
		# CDS	Trehalose-6-phosphate hydrolase	lcl|contig00186:126780-128441	yberc_1530
		# rRNA	ribosomal RNA	lcl|contig00186:128108->128665	yberc_r60
		if ( $g2start <=  $g1end && $g1ctg eq $g2ctg 
			   && $gene2 =~ /^$id:(t|r)RNA/ && $gene1 =~ /^$id:CDS/ ) {
 				push @todelete,$gene1;
 				print "RNA overlaps CDS, +1, will remove CDS: $gene1" if $self->debug;
 		}
 		if ( $g1start <=  $g2end && $g1ctg eq $g2ctg 
			   && $gene1 =~ /^$id:(t|r)RNA/ && $gene2 =~ /^$id:CDS/ ) {
 				push @todelete,$gene2;
 				print "RNA overlaps CDS, +1, will remove CDS: $gene2" if $self->debug;
 		}
		if ( $g2end >=  $g1start && $g1ctg eq $g2ctg 
			   && $gene2 =~ /^$id:(t|r)RNA/ && $gene1 =~ /^$id:CDS/ ) {
 				push @todelete,$gene1;
 				print "RNA overlaps CDS, +1, will remove CDS: $gene1" if $self->debug;
 		}
 		if ( $g1end >=  $g2start && $g1ctg eq $g2ctg 
			   && $gene1 =~ /^(t|r)$id:RNA/ && $gene2 =~ /^$id:CDS/ ) {
 				push @todelete,$gene2;
 				print "RNA overlaps CDS, +1, will remove CDS: $gene2" if $self->debug;
 		}

		if ( $g2end >=  $g1end && $g1ctg eq $g2ctg 
			   && $gene2 =~ /^$id:(t|r)RNA/ && $gene1 =~ /^$id:CDS/ ) {
 				push @todelete,$gene1;
 				print "RNA overlaps CDS, -1, will remove CDS: $gene1" if $self->debug;
 		}
 		if ( $g1end >=  $g2end && $g1ctg eq $g2ctg 
			   && $gene1 =~ /^$id:(t|r)RNA/ && $gene2 =~ /^$id:CDS/ ) {
 				push @todelete,$gene2;
 				print "RNA overlaps CDS, -1, will remove CDS: $gene2" if $self->debug;
 		}

		$gene1 = $gene2;
	}

	@todelete = unique(@todelete);
	@todelete;
}

# DiscRep_ALL:RNA_CDS_OVERLAP::6 coding regions overlap RNA features
# DiscRep_SUB:RNA_CDS_OVERLAP::1 coding regions are completely contained in RNAs
# bcere0010:CDShypothetical proteinlcl|contig00145:2461-2892bcere0010_53360
# bcere0010:rRNA23S ribosomal RNAlcl|contig00145:9-2928bcere0010_r50
sub get_from_discrp {
	my ($self,$header) = @_;
	my @lines;
	my $readflag = 0;

	open MYIN,"discrp" or die "Could not open discrp file";

	while (<MYIN>) {
		$readflag = 0 if ( /^\n/ || /^\s/ );
		push @lines,$_ if $readflag;
		$readflag = 1 if ( /DiscRep_SUB:$header/ );
		print "Found header $header in discrp\n" if ($self->debug && $readflag);
		# Older versions of tbl2asn have a slightly different 'discrp' format
		# $readflag = 1 if ( /^$header/ );
	}
	@lines;
}

sub write_tbl {
	my ($self,$tbl) = @_;

	my $dir = $self->outdir;
	my $id = $self->id;
	my $namemap = $self->namemap;
    # my $count = 1;

	`mv $dir/$id.tbl $dir/$id.tbl.orig` if ( -e "$dir/$id.tbl" );

    my $tblfh = FileHandle->new(">$dir/$id.tbl") or die ("Cannot write to file $dir/$id.tbl");

	print "Writing to $dir/$id.tbl\n" if $self->debug;

	for my $contig ( @{$tbl} ) {

		#print $tblfh ">Features " . $contig->{contigname} . " Table${count}\n";
        print $tblfh ">Features " . $contig->{contigname} . "\n";

	 FEAT:
		for my $feat ( sort sort_by_loc keys %{$contig} ) {
			next FEAT if ( $feat eq 'contigname' );
			print $tblfh $feat;
			print $tblfh $contig->{$feat};
		}
        #$count++;
	}

	1;
}

sub is_short {
	my ($self,$contigname) = @_;
	my $namemap = $self->namemap;
	my $cutoff = $self->cutoff;

	for my $name ( keys %{$namemap} ) {
		if ( $contigname =~ /$name/ && $namemap->{$name}->{len} < $cutoff ) {
			print "Skipping contig " . $contigname . ", length:" . 
			  $namemap->{$name}->{len} . "\n" if $self->debug;
			return 1;
		}
	}
	0;
}

sub sort_by_loc {
	my ($a1) = $a =~ /^[<>]?(\d+)/; 
	my ($b1) = $b =~ /^[<>]?(\d+)/; 

	$a1 <=> $b1;
}

sub delete_from_tbl {
	my ($self,$tbl,@todelete) = @_;
	my %todelete;

	# Transform array of lines into hash of names plus locations
	# Gene    yberc_r90       lcl|contig00186:128269->128665  yberc_r90
	# CDS     hypothetical protein    lcl|contig00890:c203-<1 yberc_39200
	for my $delete (@todelete) {
		$delete =~ /^\S+\s+[^|]+\|((?:contig|ctg)[\d.]+):c?(\d+)[<>]?-[<>]?(\d+)/;
		$todelete{"$1 $2 $3"}++;
	}

	# Iterate over each contig section
	for my $contig ( @{$tbl} ) {
		print "Next contig in *tbl is " . $contig->{contigname} . "\n" 
		  if $self->debug;

		for my $delete ( keys %todelete ) {
			my ($deletecontig,$loc1,$loc2) = $delete =~ /^(\S+)\s(\d+)\s(\d+)/;

			if ( $deletecontig eq $contig->{contigname} ) {
				for my $feat ( keys %{$contig} ) {
					my ($loca,$locb) = $feat =~ /^[<>]?(\d+)\s+[<>]?(\d+)/;
					if ( $loca == $loc1 && $locb == $loc2 ) {
						print "Deleting feature $delete with location \'$loc1 $loc2\'\n" 
						  if $self->debug;
						delete ${$contig}{$feat};
					}
				}
			}

		}
	}

	$tbl;
}

# Each contig, starting with ">Features", is an element in
# the array, the features and the name of the contig
# are stored in an anonymous hash, one per element.
# Also store the contig names in an array.
sub read_tbl {
	my $self = shift;
	my $tbl = $self->outdir . "/" . $self->id . ".tbl";
	my (@contignames,@tbl);

	local $/ = undef;

	open MYIN,"$tbl" or die "Cannot open file $tbl";

	my @contigs = split  /^(?=>Features.+)/m, <MYIN>; 

	for my $contig ( @contigs ) {

		my $hsh;
		$contig =~ /^>Features\s+(\S+)/;
		$hsh->{contigname} = $1;
		push @contignames,$1;

		my @features = split /^(?=[\d><]+\s+[\d><]+\s+\w+)/m, $contig;
		# Remove the header
		my $header = shift @features;
		print "Header: $header\n" if $self->debug;

		for my $feature ( @features ) {
			$feature =~ /^([\d><]+\s+[\d><]+\s+\S+)(.+)/s;
         # print Dumper $feature;
			$hsh->{$1} = $2;
		}
      push @tbl, $hsh;
	}
	$self->contigs(\@contignames);

	\@tbl;
}

sub cleanup {
	my $self = shift;
	my $id = $self->id;
	my $dir = $self->outdir;
	my $template = $self->{template};

	`rm $template.sbt` if -e "$template.sbt";
	`rm discrp.orig`   if -e "discrp.orig";
	`mv discrp $dir`   if -e "discrp";
        `mv $id.agp $dir`  if -e "$id.agp";

	for my $suffix ( qw(gbf val tbl sqn) ) {
		unlink "$dir/$id.$suffix.orig" if -e "$dir/$id.$suffix.orig";
	}

}

sub outdir {
	my ($self,$id) = @_;

	if ($id) {
		# Make output directory
		my $outdir = "$id-gbsubmit.out.d";
		`mkdir -p $outdir` if ( ! -d $outdir );
		$self->{outdir} = $outdir;
	}
	return $self->{outdir} if $self->{outdir};
}

sub readsPerBase {
	my ($self,$len,$avg) = @_;

	$self->{readsPerBase}->{$len} = $avg if ($len && $avg);

	if (defined $self->{readsPerBase}) {
		return $self->{readsPerBase};
	} else {
		return 0;
	}
}

sub edit_asn_file {
    my $self = shift;

# Put species into the title, e.g. replace "Annotation of the XXXXX YYYYY genome"
# with "Annotation of the Yersinia bercovieri ATCC 12345 genome"
    my $species = $self->organism;
    my $strain  = $self->strain;

    local $/ = undef;

    my $template = $self->template;
    my $asn      = $template . ".sbt";
    unlink $asn if ( -e $asn );

    open MYIN, $template or die "Cannot open file $template";
    my $text = <MYIN>;
    $text =~ s/<BACTERIALSPECIES>/$species/;
    $text =~ s/<BACTERIALSTRAIN>/$strain/;

    open MYOUT, ">$asn" or die "Cannot write to file $asn";
    print MYOUT $text;

}

sub make_top_comment {
    my $self = shift;

    my $comment = "Bacteria available from MARC" if ( $self->strain !~ /ATCC/ );

    $comment;
}

sub create_qual {
    my ( $self, $qualfile ) = @_;

    my $contignames = $self->contigs;
    my $outdir      = $self->outdir;
    my $id          = $self->id;
    my %qualities   = ();

    if ( -e $qualfile ) {
        my $in = Bio::SeqIO->new(
            -file   => $qualfile,
            -format => "qual"
        );

        while ( my $qual = $in->next_seq ) {
            my $primary_id = $qual->primary_id;
            $qualities{$primary_id} = $qual;
        }

        my $out = Bio::SeqIO->new(
            -file   => ">>$outdir/$id.qvl",
            -format => "qual"
        );

        # A name could be something like "contig01112.6" or "contig01112"
        for my $name ( @{$contignames} ) {

            for my $key ( keys %qualities ) {
                if ( $name =~ /$key/ ) {
                    my $obj = $qualities{$key};
                    $obj->display_id($name);
                    $out->write_seq($obj);
                }
            }
        }

    }
    else {
        print "No file $qualfile found\n" if $self->debug;
    }
}

# # ORGANISM: Homo sapiens
# # TAX_ID: 9606
# # ASSEMBLY NAME: EG1
# # ASSEMBLY DATE: 06-September-2006
# # GENOME CENTER: NCBI
# # DESCRIPTION: Example AGP specifying the assembly of chromosome Y from WGS contigs
# chrY 1 3043 1 W AADB02037551.1 1 3043 +
# spacer: 'Chr1\t<start>\t<end>\t<number>\tN\t100\tfragment\tyes'
sub create_agp {
    my ( $self, $gbk ) = @_;
    my $taxid    = $self->taxid or die "No taxid found";
    my %frames   = ( '1', '+', '-1', '-' );

    my ($date,$id,$organism,$strain) = 
    ($self->get_date,$self->id,$self->organism,$self->strain);

# hack!
# my @gbks = <*rnammer.out.gbk>;
# die "No *rnammer.out.gbk file found, can not create an *.agp file" unless ( -e $gbks[0] );

    open MYAPG, ">>$id.agp" or die "Cannot create file $id.agp";

    my $in = Bio::SeqIO->new( -file => $gbk, -format => 'genbank' );
    my $seq = $in->next_seq;
    my @contigs =
      grep { $_->primary_tag eq 'fasta_record' } $seq->get_SeqFeatures;

    my $text = "# ORGANISM: $organism
# TAX_ID: $taxid
# ASSEMBLY NAME: $id
# ASSEMBLY DATE: $date
# GENOME CENTER: NMRC
# DESCRIPTION: $organism $strain chromosome, whole genome shotgun
";

    print MYAPG $text;

    my $count      = 1;
    my $pos        = 1;
    my $spacer_len = 100;

    for my $contig (@contigs) {
        my @names = $contig->get_tag_values('name');
        my $len = ( $contig->end ) - ( $contig->start );
        print MYAPG "Chr1\t" 
          . $pos . "\t"
          . ( $len + $pos ) . "\t"
          . $count . "\tW\t"
          . $names[0] . "\t1\t"
          . ( $len + 1 ) . "\t"
          . ( $frames{ $contig->strand } ) . "\n";
        $pos = $len + $pos + 1;
        $count++;
        print MYAPG "Chr1\t" 
          . $pos . "\t"
          . ( $pos + $spacer_len - 1 ) . "\t"
          . $count
          . "\tN\t$spacer_len\tfragment\tyes\n";
        $pos += $spacer_len;
        $count++;
    }

}

# The *cmt file is a 2 column, tab-delimited file like this:
# StructuredCommentPrefix ##Genome-Assembly-Data-START##
# Assembly Method Newbler v. 2.3
# Assembly Name   Ecoli.str12345_v1.0
# Genome Coverage 16.3x
# Sequencing Technology   454 Titanium; PacBio RS
# StructuredCommentSuffix ##Genome-Assembly-Data-END##
sub create_cmt {
    my $self = shift;
    my ( $method, $name, $coverage, $tech, $id, $outdir, $readsPerBase ) = (
        $self->Assembly_Method, $self->Assembly_Name,
        $self->Genome_Coverage, $self->Sequencing_Technology,
        $self->id,              $self->outdir,
        $self->readsPerBase
    );
    $method = trim_comma($method);
    $tech   = trim_comma($tech);
    my ( $totalLen, $totalReads, $cmtfh );

    if ($readsPerBase) {
        for my $len ( keys %{$readsPerBase} ) {
            $totalLen   += $len;
            $totalReads += ( $len * $readsPerBase->{$len} );
        }
        $coverage = sprintf( "%.1f", ( $totalReads / $totalLen ) );
    }

    if ( -e "$outdir/$id.cmt" ) {
      $cmtfh = FileHandle->new(">>$outdir/$id.cmt")
        or die("Cannot open file $outdir/$id.cmt for appending");
    } else {
      $cmtfh = FileHandle->new(">$outdir/$id.cmt")
        or die("Cannot open file $outdir/$id.cmt for writing");      
    }

    my $txt = "StructuredCommentPrefix\t" . '##Genome-Assembly-Data-START##' . "\n" .
              "Assembly Method\t$method\n";
    $txt .=   "Assembly Name\t$name\n" if $name;
    $txt .=   "Genome Coverage\t${coverage}x\n" .
              "Sequencing Technology\t$tech\n" .
              "StructuredCommentSuffix\t" . '##Genome-Assembly-Data-END##';

    print $cmtfh $txt;
    1;
}

sub get_date {
	my $self = shift;
	return time2str("%d-%B-%Y",time);
}

sub trim_comma {
    my $str = shift;
    $str =~ s/('|")$//;
    $str =~ s/^('|")//;
    $str;
}

sub edit_definition {
    my ( $self, $desc ) = @_;

 # Edit DEFINITION line, input looks something like:
 # [organism=Yersinia bercovieri] [strain=ATCC_43970] [gcode=11] [date=7-3-2008]
    my ($organism) = $desc =~ /\[organism=([^]]+)/;
    my ($strain)   = $desc =~ /\[strain=([^]]+)/;
    my ($gcode)    = $desc =~ /\[gcode=([^]]+)/;
    my $definition = '';

    if ( $strain eq undef ) {
        $strain = '';
    }
    else {
        $strain =~ s/_/ /g;
    }

    $self->organism($organism);
    $self->strain($strain);

    if ( $strain =~ /ATCC/ ) {
        my ($id) = $strain =~ /ATCC\s+(.+)/;
        $definition =
"[organism=$organism][strain=$strain][culture_collection=ATCC:$id][gcode=$gcode][tech=wgs]";
    }
    else {
        $definition =
          "[organism=$organism][strain=$strain][gcode=$gcode][tech=wgs]";
    }
    $definition;
}

# Get names, see if any are duplicated - /name="contig01098"
#     fasta_record    1..529
#                     /name="contig00231"
#                     /note="untiled"
#                     /note="Calculated sequence coverage = 5.36 reads per base"
sub make_namemap {
    my ( $self, $file ) = @_;
    my $namemap;

    my $in = Bio::SeqIO->new( -file => $file, -format => 'genbank' );
    my $seq = $in->next_seq;
    for my $feat ( $seq->get_SeqFeatures ) {
        if ( $feat->primary_tag eq 'fasta_record' ) {
            my @names = $feat->get_tag_values('name');
            $namemap->{ $names[0] }->{num}++;

            # $namemap{$1}++ if ( /\/name="([^"]+)/ );
            $namemap->{ $names[0] }->{len} = ( $feat->end ) - ( $feat->start );
        }
    }

    $self->namemap($namemap);
}

sub number_the_duplicate {
	my ($self,$name) = @_;

	if ( ! defined $self->{duplicates}->{$name} ) {
		$self->{duplicates}->{$name} = 1;
	} else {
		$self->{duplicates}->{$name}++;
	}
	($name . "." . $self->{duplicates}->{$name});
}

sub is_duplicate_name {
	my ($self,$name) = @_;

	my $namemap = $self->{namemap};
	if ( $namemap->{$name}->{num} > 1 ) {
		return 1;
	} else {
		return 0;
	}
}

sub _initialize {
	my $self = shift;

	while ( @_ ) {
		( my $key = shift ) =~ s/^-//;
		$self->{$key} = shift;
	}
}

# $property = sub {
#     my ( $self, $name, $value ) = @_;
#
#     if ($value) {
#         $self->{$name} = $value;
#         return $value;
#     }
#     else {
#         return $self->{$name};
#     }
# };

sub lastBase {
	my ($self,$base) = @_;
    $self->{'lastbase'} = $base if defined $base;
    return $self->{'lastbase'};
}

sub accession_prefix {
	my ($self,$base) = @_;
    $self->{'accession_prefix'} = $base if defined $base;
    return $self->{'accession_prefix'};
}

sub executable {
	my ($self,$base) = @_;
    $self->{'executable'} = $base if defined $base;
    return $self->{'executable'};
}

sub namemap {
	my ($self,$base) = @_;
    $self->{'namemap'} = $base if defined $base;
    return $self->{'namemap'};
}

sub cutoff {
	my ($self,$base) = @_;
    $self->{'cutoff'} = $base if defined $base;
    return $self->{'cutoff'};
}

sub debug {
	my ($self,$base) = @_;
    $self->{'debug'} = $base if defined $base;
    return $self->{'debug'};
}

sub template {
	my ($self,$base) = @_;
    $self->{'template'} = $base if defined $base;
    return $self->{'template'};
}

sub taxid {
	my ($self,$base) = @_;
    $self->{'taxid'} = $base if defined $base;
    return $self->{'taxid'};
}

sub organism {
	my ($self,$base) = @_;
    $self->{'organism'} = $base if defined $base;
    return $self->{'organism'};
}

sub id {
	my ($self,$base) = @_;
    $self->{'id'} = $base if defined $base;
    return $self->{'id'};
}

sub strain {
	my ($self,$base) = @_;
    $self->{'strain'} = $base if defined $base;
    return $self->{'strain'};
}

sub contigs {
	my ($self,$base) = @_;
    $self->{'contigs'} = $base if defined $base;
    return $self->{'contigs'};
}

sub host {
    my ($self,$base) = @_;
    $self->{'host'} = $base if defined $base;
    return $self->{'host'};
}
    
sub country {
    my ($self,$base) = @_;
    $self->{'country'} = $base if defined $base;
    return $self->{'country'};
}

# 20-May-2012
sub collection_date {
    my ($self,$base) = @_;
    # @mon = qw(Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec);
    # printf "%4d-%s-%02d\n", $d[5]+1900, $mon[$d[4]], $d[3];
    $self->{'collection_date'} = $base if defined $base;
    return $self->{'collection_date'};
}

sub isolation_source {
    my ($self,$base) = @_;
    $self->{'isolation_source'} = $base if defined $base;
    return $self->{'isolation_source'};
}
    
sub submission_note {
    my ($self,$base) = @_;
    $self->{'submission_note'} = $base if defined $base;
    return $self->{'submission_note'};
}

sub gcode {
    my ($self,$base) = @_;
    $self->{'gcode'} = $base if defined $base;
    return $self->{'gcode'};
}

sub Assembly_Name {
    my ($self,$base) = @_;
    $self->{'Assembly_Name'} = $base if defined $base;
    return $self->{'Assembly_Name'} if defined $self->{'Assembly_Name'};
    '';
}

sub Genome_Coverage {
    my ($self,$base) = @_;
    $self->{'Genome_Coverage'} = $base if defined $base;
    return $self->{'Genome_Coverage'};
}

sub Sequencing_Technology {
    my ($self,$base) = @_;
    $self->{'Sequencing_Technology'} = $base if defined $base;
    return $self->{'Sequencing_Technology'};
}

sub Assembly_Method {
    my ($self,$base) = @_;
    $self->{'Assembly_Method'} = $base if defined $base;
    return $self->{'Assembly_Method'};
}

sub trim {
	my $str = shift;
	$str =~ s/^\s+//;
	$str =~ s/\s+$//;
	$str =~ s/\s+/ /g;
	$str;
}

sub unique {
	my @arr = @_;
	my %hsh;

	for my $x (@arr) {
		$hsh{$x}++;
	}

	(keys %hsh);
}

1;

__END__


