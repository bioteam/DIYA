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
the application that creates ASN.1 for Genbank. Also does any QC on the 
CDS and gene features as required by NCBI.

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


sub new {
	my ($module,@args) = @_;
	my $class = ref($module) || $module;
	my $self = {};

	bless($self,$class);

	# defaults
	$self->template('template');
	$self->executable('/usr/local/bin/tbl2asn');

	$self->_initialize(@args);

	$self;
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
	my ($self,$feat) = @_;

	if ( $feat->primary_tag eq 'CDS' ) {

		# edit product tag and a note tag
		if ( $feat->has_tag('product') ) {

			my @product = $feat->remove_tag("product");

			# if 'note' and 'product' are duplicated
			if ( $feat->has_tag('note') ) {
				my @notes = $feat->remove_tag('note');
				for my $note (@notes) {
					$note =~ s/\s+$//;
					my $str = "similar to " . $product[0];
					$feat->add_tag_value('note',$note) unless
					  ( $note eq $str );
				}
			}

			$product[0] =~ s/^\s*//;

			# If the product comes from UniRef, something like:
			# "UPF0076 protein yjgF n=146 Tax=Bacteria RepID=YJGF_ECOL6" or
			# "Putative Orf27; P2 LysB homolog; control of lysis [Ente. n=2 Tax=Yersinia RepID=Q66BL7_YERPS"
			if ( $product[0] =~ /(.+?)\s+n=\d+\s+Tax=(.+)/ ) {

				my ($newProduct,$species) = ($1,$2);
				$species =~ s/RepID/UniRef RepID/;

				my $note = "similar to $newProduct of $species";
				$note =~ s/\sn=\d+\s/ /;
				# 'protein protein'
				$note =~ s/protein protein/protein/;
				$note =~ s/similar to Similar/similar/i;

				# Remove COG or FOG ids found in the UniRef headers, e.g.:
				# FOG: TPR repeat n=1 Tax=Vibrio vulnificus RepID=Q8DF47_VIBVU
				# COG0784: FOG: CheY-like receiver n=1 Tax=Bacillus anthracis RepID=UPI00
				$newProduct =~ s/^(FOG:\s+|COG\d+:\s+FOG:\s+|COG\d+:\s*)//;

				$feat->add_tag_value('note',$note);

				# TIGR ids
				if ( $newProduct =~ /^UPF\d+\s+(zinc-binding protein|ATP-binding protein)/ ) {
					$newProduct = $1;
				}

				# Add back the product but do not use locus tags from other 
				# genomes, e.g. product="hypothetical protein YPO0973
				if ( $newProduct =~ /^hypothetical protein\s+.*/i ) {
					$product[0] = 'hypothetical protein';
				} else {
					$product[0] = $newProduct;
				}

				$feat->remove_tag('score') if $feat->has_tag('score');

				# 'protein protein'
				$product[0] =~ s/protein protein/protein/;

			}

			# if product name comes from match to COG, not UniRef
			if ( $product[0] =~ /^(COG\d+):\s+(.+)/ ) {
				$product[0] = $2;
				my $note = "similar to $1";
				$feat->add_tag_value('note',$note) if (! $1 eq 'family' );

				#$product[0] =~ s/^Predicted/hypothetical protein/;
			}

			# if product looks something like:
			# "GpL [Enterobacteria phage P2] ..." then correct
			if ( $product[0] =~ /(\w+)\s+\[(Enterobacteria[^]]+)\]/ ) {
				$product[0] = "$2 $1";
			}

			# if the product looks something like:
			# "ABC protein [Yersinia pestis CO92]." then correct
			if ( $product[0] =~ /(.+)\s+(\[.+?\])\.$/ ) {

				my ($newProduct,$species) = ($1,$2);

				# Add back the product but do not use locus tags from other 
				# genomes, e.g. product="hypothetical protein YPO0973
				if ( $newProduct =~ /^hypothetical protein\s+.*/i ) {
					$product[0] = 'hypothetical protein';
				} else {
					$product[0] = $newProduct;
				}

				my @score = $feat->remove_tag('score') if $feat->has_tag('score');

				if ( $score[0] ) {
					my $note = "similar to $newProduct $species";
					$feat->add_tag_value('note',$note);

					#$note = "BLAST score=" . $score[0];
					#$feat->add_tag_value('note',$note);
				}
			}

			# no species names allowed in product/
			if ( $product[0] =~ /^Similar to\s+(.+?)\s+of\s+Escherichia\s+coli/ ) {
				my $note = $product[0];
				$product[0] = $1;
				$product[0] =~ s/^(putative|probable)\s+//;
				$feat->add_tag_value('note',$note);
			}

			if ( $product[0] =~ /^(similar|similarities) (to|with) (unknown|putative|probable|C-terminal)/i ) {
				$feat->add_tag_value('note', $product[0]);
				$product[0] = 'hypothetical protein';
			}

			$product[0] = correct_spelling($product[0]);

			$product[0] = add_trailing($product[0]);

			$product[0] = remove_trailing($product[0]);

			$product[0] = make_singular($product[0]);

			$product[0] = remove_banned($product[0]);

			$product[0] = 'hypothetical protein' if ( is_hypothetical($product[0]) );

			$product[0] = remove_similar($product[0]);

			$product[0] = remove_loci($product[0]);

			# Too long, too many ()'s
			if ( $product[0] =~ /^([^()]+ \([^)]+\)) \([^)]+\)/ ) {
				$product[0] = $1;
			}

			# remove '...-like'
			if ( $product[0] =~ /-like$/i ) {
				$feat->add_tag_value('note', $product[0]);
				$product[0] = 'hypothetical protein';
			}

			# remove clauses with semicolons from 'product', add to a 'note'
			if ( $product[0] =~ /^([^;]+);\s*(.+)/ ) {
				$product[0] = $1;
				$feat->add_tag_value('note', $2);
			}

			# finally add it back
			$feat->add_tag_value('product', $product[0]);

		} else {

			# if there is no 'product' then it's a 'hypothetical...'
			$feat->add_tag_value('product','hypothetical protein');

		}

		# change locus_tag to protein_id
		if ( $feat->has_tag('locus_tag') ) {
			my @loci = $feat->remove_tag('locus_tag');
			my $protein_id = 'gnl|nmrcread|' . $loci[0];
			$feat->add_tag_value('protein_id', $protein_id);
		}

		# change xref tags to note tags
		if ( $feat->has_tag('xref') ) {
			my @xrefs = $feat->remove_tag('xref');

			for my $xref ( @xrefs ) {
				my $note = "similar to $xref" ;
				$feat->add_tag_value('note',$note) if ( ! $xref eq 'family' );
			}
		}

		# features with 'rps_gi' tag need to be rearranged
		if ( $feat->has_tag('rps_gi') ) {
			my @ids = $feat->remove_tag('rps_gi') ;

			for my $id ( @ids ) {
				my $note = "similar to motif $id";
				$feat->add_tag_value('note',$note) unless ( $id eq 'family' );
			}

			# these matches have no other information
			if ( $ids[0] =~ /^(pfam|COG|cd)/ ) {
            $feat->remove_tag('product');
            $feat->add_tag_value('product','hypothetical protein');
			}

			# change cluster, score, and group tags to a note tag if there's data
			my @clusters = $feat->remove_tag('cluster') if $feat->has_tag('cluster');
			my @groups = $feat->remove_tag('group') if $feat->has_tag('group');
			my @scores = $feat->remove_tag('score') if $feat->has_tag('score');

			if ( $clusters[0] && $groups[0] && $scores[0] ) {
				my $note = "similar to CDD " . $clusters[0] . ", group " . $groups[0];
				$feat->add_tag_value('note',$note);

				#$note = "RPS BLAST score=" . $scores[0];
				#$feat->add_tag_value('note',$note);
			}
		}

                # remove all scores
                if ( $feat->has_tag('score') || $feat->has_tag('Score') ) {
                   my @scores = $feat->remove_tag('score') ;
                }
 
                for my $tag ($feat->get_all_tags) {                    
                   for my $value ($feat->get_tag_values($tag)) {
		       $feat->remove_tag($tag) if ( $value =~ /^\s*score\s*[=:]/i );
                   }          
                }       

		return $feat;
	}

	elsif ( $feat->primary_tag eq 'gene' ) {

		return $feat;
	}

	elsif ( $feat->primary_tag eq 'tRNA' ) {

		# add inference tag
		$feat->add_tag_value('inference','profile:tRNAscan-SE:1.23') 
		  if ( ! $feat->has_tag("inference") );

		# Edit note tag
		if ( $feat->has_tag("note") ) {
			my @note = $feat->remove_tag("note");
			$note[0] =~ s/score=/tRNAscan-SE score=/;
			$feat->add_tag_value('note', $note[0]);
		}

		# Change Codon tag to a note tag
		if ( $feat->has_tag("Codon") ) {
			my @codon = $feat->remove_tag("Codon");
			$feat->add_tag_value('note',"Anticodon is $codon[0]");
		}

		# remove AminoAcid tag
		$feat->remove_tag("AminoAcid") if $feat->has_tag('AminoAcid');

		# remove 'ID', which will be added back as 'product'
		if ( $feat->has_tag("ID") ) {
			my @ID = $feat->remove_tag("ID");
			$ID[0] =~ s/:/-/;
			$feat->add_tag_value('product',$ID[0]);
		}

		# add a 'gene' feature for the tRNA
		my $genefeat = Bio::SeqFeature::Generic->new(-primary_tag => 'gene');
		$genefeat->location($feat->location);

		($genefeat->end) > ($genefeat->start) ? $genefeat->strand(1) :
		  $genefeat->strand(-1);

		my @loci = $feat->remove_tag('locus_tag') if $feat->has_tag('locus_tag');
		$genefeat->add_tag_value('locus_tag', $loci[0]);

		return ($feat,$genefeat);
	}

	elsif ( $feat->primary_tag eq 'rRNA' ) {

		# delete note tag with "frame=." value
		if ( $feat->has_tag('note') ) {
			my @notes = $feat->remove_tag("note") ;

			for my $note (@notes) {
				$note =~ s/score=/RNAMMER score=/;
				$feat->add_tag_value('note',$note) if ( $note !~ /frame=\./ );
			}
		}

		# add a 'gene' feature for the rRNA
		my $genefeat = Bio::SeqFeature::Generic->new(-primary_tag => 'gene');
		$genefeat->location($feat->location);

		($genefeat->end) > ($genefeat->start) ? $genefeat->strand(1) :
		  $genefeat->strand(-1);

		my @loci = $feat->remove_tag('locus_tag') if $feat->has_tag('locus_tag');

		$genefeat->add_tag_value('locus_tag', $loci[0]);

		return ($feat,$genefeat);
	}

}

sub correct_spelling {
	my $product = shift;

	$product =~ s/Catabolite protein activator/catabolite protein activator/;
	$product =~ s/Accessory protein regulator ([A-Z])/accessory protein regulator $1/;
	$product =~ s/Penicillin V acylase\. Cysteine peptidase\./penicillin V acylase; cysteine peptidase;/;
	$product =~ s/FtsH-2 peptidase. Metallo peptidase. MEROPS family M41,/FtsH-2 peptidase; Metallo peptidase; MEROPS family M41;/;
	$product =~ s/Camelysin\./Camelysin;/i;
	$product =~ s/Metallo peptidase\./Metallo peptidase;/i;
	$product =~ s/D-Ala-D-Ala carboxypeptidase A\./D-Ala-D-Ala carboxypeptidase A;/i;
	$product =~ s/Serine peptidase\./Serine peptidase;/i;
	$product =~ s/carboxypeptidase DacF\./carboxypeptidase DacF;/i;
	$product =~ s/acetoincleaving/acetoin cleaving/i;
	$product =~ s/dyhydrogenase/dehydrogenase/i;
	$product =~ s/carrierprotein/carrier protein/i;
	$product =~ s/proteinl/protein/ig;
	$product =~ s/POLY\(A\) POLYMERASE \/ TRNA NUCLEOTIDYLTRANSFERASE/tRNA nucleotidyltransferase/i;
	$product =~ s/MOSQUITOCIDAL TOXIN PROTEIN/mosquitocidal toxin protein/;
	$product =~ s/MONO-ADP-RIBOSYLTRANSFERASE C3/mono-ADP-ribosyltransferase C3/;
	$product =~ s/hyothetical/hypothetical/ig;
	$product =~ s/biosynthsis/biosynthesis/ig;
	$product =~ s/AMIDOHYDROLASE/amidohydrolase/;
	$product =~ s/ACETOIN TRANSPORT PERMEASE PROTEIN/acetoin transport permease protein/;
	$product =~ s/(trasporter|transoprted|tranporter)/transporter/gi;
	$product =~ s/3-OXOADIPATE ENOL-LACTONASE/3-oxoadipate enol-lactonase/;
	$product =~ s/4''''-phosphopantetheinyl/phosphopantetheinyl/;
	$product =~ s/ \(And other\) / /i;
	$product =~ s/sensor\(S\)/sensors/ig;
	$product =~ s/TRANSPOSON/Transposon/g;
	$product =~ s/INVERTASE/Invertase/g;
   $product =~ s/\bsopre\b/spore/ig;
   $product =~ s/proteintic/protein/ig;
   $product =~ s/NAD\(P\)H/NADPH/ig;
   $product =~ s/spaning/spanning/ig;
	$product =~ s/proetin/protein/ig;
	$product =~ s/fibre/fiber/ig;
	$product =~ s/Hypotethical/hypothetical/ig;
	$product =~ s/reguatory/regulatory/ig;
	$product =~ s/fibre/fiber/ig;
	$product =~ s/Uncharacterised/Uncharacterized/ig;
	$product =~ s/putaive/putative/ig;
	$product =~ s/haemin/hemin/ig;
	$product =~ s/Haemolytic/Hemolytic/ig;
	$product =~ s/unknow /unknown /ig;
	$product =~ s/haemagglutinin/hemagglutinin/ig;
	$product =~ s/PERIPLASMIC PROTEIN/periplasmic protein/g;
	$product =~ s/MITOCHONDRIAL TRANSPORTER ATM1 \(Atm1\)/Mitochondrial transporter Atm1/;
	$product =~ s/Protein C\. Serine peptidase\./Protein C serine peptidase/;
	$product =~ s/unkown/unknown/ig;
	$product =~ s/addtional/additional/ig;

	$product;
}

sub remove_loci {
	my $product = shift;

	# no loci in product names, e.g 'hydrolase LC_123'
	return $1 if ( 
			$product =~ /^(Uncharacterized adenine-specific methylase)\s+[\d{1,}a-z{1,}\/-_]+$/i ||
         $product =~ /^(hydrolase)\s+[\d{1,}a-z{1,}\/-_]+$/i ||
			$product =~ /^(metalloprotease)\s+[\d{1,}a-z{1,}\/-_]+$/i ||
			$product =~ /^(Phosphodiesterase),?\s+[\d{1,}a-z{1,}\/-_]+$/i ||
			$product =~ /^(reductase)\s+[\d{1,}a-z{1,}\/-_]+$/i ||
			$product =~ /^(transport protein)\s+[\d{1,}a-z{1,}\/-_]+$/i ||
         $product =~ /^(phosphotransferase)\s+[\d{1,}a-z{1,}\/-_]+$/i ||
         $product =~ /^(Uncharacterised conserved protein)\s+[\d{1,}a-z{1,}\/-_]+$/i ||
         $product =~ /^(Uncharacterized MFS-type transporter)\s+[\d{1,}a-z{1,}\/-_]+$/i ||
         $product =~ /^(Uncharacterized transporter)\s+[\d{1,}a-z{1,}\/-_]+$/i ||
         $product =~ /^(Pirin-like protein)\s+[\d{1,}a-z{1,}\/-_]+$/i ||
         $product =~ /^(Uncharacterized protease)\s+[\d{1,}a-z{1,}\/-_]+$/i ||
         $product =~ /^(HAD-superfamily hydrolase, subfamily IB,)\s+[\d{1,}a-z{1,}\/-_]+$/i ||
         $product =~ /^(Uncharacterized mscS family protein)\s+[\d{1,}a-z{1,}\/-_]+$/i ||
         $product =~ /^(Shikimate 5-dehydrogenase-like protein)\s+[\d{1,}a-z{1,}\/-_]+$/i ||
         $product =~ /^(Peptidase T-like protein)\s+[\d{1,}a-z{1,}\/-_]+$/i ||
         $product =~ /^(Uncharacterized ABC transporter ATP-binding protein)\s+[\d{1,}a-z{1,}\/-_]+$/i ||
         $product =~ /^(Uncharacterized acyl-CoA thioester hydrolase)\s+[\d{1,}a-z{1,}\/-_]+$/i
        );

	$product;
}

sub add_trailing {
	my $product = shift;
	
	$product .= "-containing protein" if ( $product =~ /^\s*\w+\s+repeat\s*$/i );

   $product =~ s/transcriptional regulator, PadR-like\s*$/transcriptional regulator, PadR-like protein/;	
	$product =~ s/DinB family superfamily\s*$/DinB superfamily protein/;
	$product =~ s/YfhF-\s*$/YfhF-like/;
	$product =~ s/WbqC-\s*$/WbqC-like/;
	$product =~ s/Zinc finger, CHC2-\s*$/Zinc finger, CHC2-like/;
   $product =~ s/Tetratricopeptide TPR_4\s*$/Tetratricopeptide TPR_4 containing protein/;
   $product =~ s/sensor histidine kinase domain\s*$/sensor histidine kinase domain containing protein/;
   $product =~ s/SWIM zinc finger\s*$/SWIM zinc finger containing protein/;
   $product =~ s/Tetratricopeptide \(TPR\) repeat\s*$/Tetratricopeptide \(TPR\) repeat containing protein/;
	$product =~ s/transporter related$/transporter related protein/i;
	$product =~ s/^(\w+) binding domain$/$1 binding domain containing protein/;
	$product =~ s/(\S+-like)$/$1 protein/;

	$product;
}

sub remove_trailing {
	my $product = shift;

	# remove trailing non-text
	$product =~ s/\s+$//;
	$product =~ s/(_|,|\?)$//;

	# Remove meaningless trailing comments
	$product =~ s/\s+and other pencillin-binding proteins?//i;
	$product =~ s/\s+\(Insoluble fraction\)$//i;
   $product =~ s/\s+\(amino terminus\)$//i;
	$product =~ s/,?\s+(SA2311|A118|AAA_5|MJ\d+|YLL\d+|SAB\d+[a-z]+|alr\d+)$//i;
	$product =~ s/\/Dioxygenas$//i;
	$product =~ s/s+domain 1 \(GGDEF\)//i;
	$product =~ s/\s+C X region$//i;
	$product =~ s/\s+\(subunit\)//i;
	$product =~ s/\s+firmicutes//i;
	$product =~ s/\(Dihydro>//i;
   $product =~ s/\s+in\S+\s+\S+region$//i;
	$product =~ s/\s+catalytic region//i;
	$product =~ s/^Acyl transferase region$/Acyl transferase/i;
	$product =~ s/permease-associated region$/permease domain protein/i;
	$product =~ s/\s+and\s+(dehydrogenase|aminotransferase)$//i;
   $product =~ s/s+\(Glutamate-aspartate carrierprotein\)$//i;
   $product =~ s/\s+\(Replication protein ori\d+\)$//i;
   $product =~ s/\s+\(Replication protein\)$//;
   $product =~ s/ in Marinococcus halophilus$//i;
   $product =~ s/\s+clpC\/mecB$//i;
   $product =~ s/\s+CA_[A-Z]+\d+$//i;
	$product =~ s/,?\s+fhuD$//i;
   $product =~ s/\s+BCE_\d+$//i;
   $product =~ s/\s+BLi\d+\/BL\d+$//i;
   $product =~ s/\s+(BT|LMOf)\d+_\d+$//i;
	$product =~ s/family protein$//i;
   $product =~ s/ HD sub domain$//i;
   $product =~ s/\s*(SA|BH|VC|NMB|HI|SH)\d+$//i;
   $product =~ s/\s+pXO\d+-\d+\/BXB\d+\/GBAA_pXO\d+_\d+$//i;
   $product =~ s/\s+BA_\d+\/GBAA\d+\/BAS\d+$//i;
   $product =~ s/\s*\(Ans operon repressor\s*protein\)$//i;
   $product =~ s/\s*\(Divided with OB2865 and OB2866\)$//i;
   $product =~ s/\s*\(N-terminal domain to N-Acetylmuramoyl-L-alanine amidase\)$//i;
   $product =~ s/\s+family\sfamily$//i;
	$product =~ s/\s+related$//i;
   $product =~ s/,? family$//;
   $product =~ s/\s*of V. anguillarum \(YclQ protein\)$//i;
	$product =~ s/\s*\(Putative endopeptidase inhibitor\)$//;
	$product =~ s/\s*\(Putative arconitate hydratase\)$//;
	$product =~ s/\s*\(Two component system response regulatory protein\)$//;
   $product =~ s/\s*\(Hypothetical protein\)$//;
	$product =~ s/\s*\(Outer membrane usher protein\)$//;
	$product =~ s/\s*\(eIF-2Bgamma.eIF-2Bepsilon\)$//;
	$product =~ s/\s*\(some contain LysM.invasin domains\)$//;
	$product =~ s/, catalytic domain:D-isomer specific 2-hydroxyacid dehydrogenase, NAD binding domain$//i;
	$product =~ s/,? truncat(ion|ed)$//i;
	$product =~ s/,? interruption$//i;
	$product =~ s/,? YHCS B.subtilis ortholog$//i;
	$product =~ s/,? \([A-Z\d-]+\)$//i;
	$product =~ s/, similar to SW:\w+$//i;
	$product =~ s/,? hexapeptide repeat$//i;
	$product =~ s/,? (C|N)-termin(us|al)$//i;
	$product =~ s/,? (N|C)-terminal (domain|region)$//i;
	$product =~ s/,? (SHL|HD) domains?$//i;
   $product =~ s/ HD sub domain$//i;
	$product =~ s/,?\s+N-region$//i;
	$product =~ s/ and inactivated derivatives$//i;
   $product =~ s/\s+and\s+orf\d+$//i;
	$product =~ s/\s+and$//;
   $product =~ s/,? (mitochondrial|chloroplastic)//;
   $product =~ s/\s*BCE?_\d+$//;
	$product =~ s/\(Dihydrolipoamide acetyltransferase component of acetoin cleavingsystem\)$//i;
   $product =~ s/\s*\(Acetoin dehydrogenase E2 component\)$//i;
   $product =~ s/\s*, small chain VC\d+$//;
   $product =~ s/\s*, GNAT family family protein$//;
   $product =~ s/\s*BC_\d+$//;
   $product =~ s/\s*, periplsmic ligand-binding protein$//;
   $product =~ s/\s*BA_\d+\/GBAA\d+\/BAS\d+$//;
   $product =~ s/\s*(\[NAD\(P\)+\]|\[Mn\]|\[NAD\(P\)H\]|\[NADPH?\]|\[ATP\]|\[Cu-Zn\]|\[ATP hydrolyzing\]|\[isomerizing\]|\[carboxylating\]|\[glutamine-hydrolyzing\]|\[decarboxylating\]|\[a?symmetrical\])$//i;
	$product =~ s/\/isomerase:Polysaccharide biosynthesis protein CapD:dTDP-4-dehydrorhamnose reductase:Nucleotide sugar epimerase$//i;
	$product =~ s/\s+\(Diaminohydroxyphosphoribosylaminopyrimidine deaminase \(Riboflavin-specific deaminase\) and 5-amino-6-\(5-phosphoribosylamino\)uracil reductase\)$//i;
	$product =~ s/\s+\(Probable\), IS891\/IS1136\/IS1341:Transposase, IS605 OrfB$//i;
	$product =~ s/\s+and inactivated derivatives-like protein$//i;

	$product;
}

sub remove_similar {
	my $product = shift;

	# remove phrases like 'similar to...'
	return $1 if $product =~ /^Similar to (acetyltransferase|exochitinase|restin isoform b|permease protein of ABC transport system|protein gp49 from prophage N15|fimbrial subunit type 1|bacteriophage integrase|vgrG protein|base plate protein gp25 of Bacteriophage|bacteriophage tail fiber assembly protein|protein V)/i;
	return $1 if $product =~ /^Related to (galactoside O-acetyltransferase)/;
	return $1 if $product =~ /^Strongly similar to (S-adenosylmethionine synthetase)/;
	return $1 if $product =~ /^Exhibited considerable similarity to a small (polypeptide \(RepA\))/;
	return $1 if $product =~ /^Similarities with (transcription regulator LysR family|alpha replication protein of prophage|Photorhabdus cytotoxin|tail fiber protein)/i;

	$product;
}

sub make_singular {
	my $product = shift;

	$product =~ s/GTPases/GTPase/ig;
	$product =~ s/variants/variant/ig;
   $product =~ s/synthetases/synthetase/ig;
	$product =~ s/glucosidases/glucosidase/ig;
	$product =~ s/thioredoxins/thioredoxin/ig;
   $product =~ s/asparaginases/asparaginase/ig;
	$product =~ s/acetylases/acetylase/ig;
	$product =~ s/enzymes/enzyme/ig;
	$product =~ s/Flavodoxins/Flavodoxin/ig;
	$product =~ s/toxins/toxin/ig;
	$product =~ s/Permeases/Permease/ig;
	$product =~ s/components/component/ig;
	$product =~ s/proteins/protein/ig;
	$product =~ s/systems/system/ig;
	$product =~ s/regulators/regulator/ig;
	$product =~ s/phosphatases/phosphatase/ig;
	$product =~ s/determinants/determinant/ig;
	$product =~ s/recombinases/recombinase/ig;
	$product =~ s/transferases/transferase/ig;
	$product =~ s/ATPases/ATPase/ig;
	$product =~ s/exporters/exporter/ig;
	$product =~ s/hydrolases/hydrolase/ig;
	$product =~ s/reductases/reductase/ig;
	$product =~ s/cytochromes/cytochrome/ig;
	$product =~ s/proteases/protease/ig;
	$product =~ s/kinases/kinase/ig;
	$product =~ s/transporters/transporter/ig;
	$product =~ s/oxidases/oxidase/ig;
	$product =~ s/helicases/helicase/ig;
	$product =~ s/synthases/synthase/ig;
	$product =~ s/peptidases/peptidase/ig;
	$product =~ s/Dehydrogenases/Dehydrogenase/ig;
	$product =~ s/lyase.+?and.+?lyases/lyase/is;

	$product;
}

sub remove_banned {
	my $product = shift;

	$product =~ s/\s+related\s+/ /ig;
	$product =~ s/\s+homologs?\s*/ /ig;
	$product =~ s/\s+\(partial\)//ig;
	$product =~ s/,?\s*putative$//i;
	$product =~ s/^(Probable|Possible|Predicted|Putative)\s+//i;
	$product =~ s/\s+\(Fragment\)\s?/ /i;
	$product =~ s/\bgene\b/\bprotein\b/;
	$product =~ s/^PREDICTED:\s*//i;
	$product =~ s/^(B.thurinienis|Salmonella)\s+//;
	$product =~ s/^Similar to\s+//i;
	$product =~ s/^Truncated\s+//i;

	$product;
}

sub is_hypothetical {
	my $product = shift;

	return 1 if ( ! $product || $product =~ /^\s+$/ );

	return 1 if (
      $product =~ /polypeptide \([a-z]+\) encoded by plasmid/i ||
      $product =~ /^Possible peptide antibiotic/ ||
      $product =~ /^PugilistDominant/i ||
      $product =~ /mutants block sporulation after engulfment/i ||
      $product =~ /^Homo sapiens/ ||
      $product =~ /^Similar to ORF13 of enterococcus faecalis/ ||
      $product =~ /^ORF13 of enterococcus faecalis TRANSPOSON TN916/i ||
      $product =~ /^CheY-homologous receiver domain / ||
      $product =~ /^Plasmid pPOD2000 Rep, RapAB, RapA, ParA, ParB, and ParC/ ||
      $product =~ /^Plasmid pRiA4b ORF-3/ ||
      $product =~ /DEHA0C09658g Debaryomyces hansenii/ ||
		$product =~ /Bacillus cereus group-specific protein/ ||
		$product =~ /^UPF\d+.+?protein.+?\S+$/i ||
		$product =~ /^BC\d+\w+ protein/i ||
  	   $product =~ /N terminal region of phage-related contractile tail sheath protein/i ||
      $product =~ /chromosome\s+\d+open\s+reading\s+frame\s+\d+/i ||
  	   $product =~ /complete genome/i ||
      $product =~ /Genome sequencing data, contig C\d+/i ||
      $product =~ /chromosome \d+ open reading frame \d+/i ||
      $product =~ /complete nucleotide sequence/i ||
      $product =~ /DNA, complete sequence/i ||
  	   $product =~ /^Genomic DNA/i ||
  	   $product =~ /whole genome shotgun sequence/i ||
  	   $product =~ /Gene, complete cds/i ||
      $product =~ /^Orf\s+\d+[A-Z]+$/i ||
		$product =~ /^Nucleoside recognition/i ||
      $product =~ /^Required for plasmid stability$/  ||
      $product =~ /^Possible Bacterial Ig-like domain/i  ||
      $product =~ /^Alpha-peptide$/i  ||
      $product =~ /LPXTG-motif cell wall anchor domain$/ ||
      $product =~ /^Amino acid tranporter$/i ||
      $product =~ /^Ankyrin repeats containing protein$/i ||
      $product =~ /^Biotin\/lipoyl attachment$/i ||
		$product =~ /^Antigen$/i ||
      $product =~ /^Restriction modification system DNA specificity domain/i ||
		$product =~ /^ABC \(ATP-binding cassette\) transporter nucleotide-binding domain$/i ||
      $product =~ /^(SAF|NACHT|Resolvase|FRG|C1q)\s+domain$/i  ||
	   $product =~ /^Contains cell adhesion domain$/i ||
		$product =~ /^Gene, IS3-like element$/i ||
      $product =~ /^Micrococcal nuclease-like protein$/ ||
      $product =~ /^modification methylase OB\d+$/i ||
      $product =~ /^Micrococcal nuclease$/i ||
      $product =~ /^leucine-rich repeat-containing protein DDB\d+$/ ||
      $product =~ /^Possible sensor histidine kinase domain/i  ||
      $product =~ /^PIN .PilT N terminus. domain$/ ||
      $product =~ /^Divergent AAA region$/i  ||
      $product =~ /^Conserved repeat domain protein$/i  ||
      $product =~ /^Collagen triple helix repeat$/i  ||
      $product =~ /^Parallel beta-helix repeat$/i  ||
      $product =~ /^PBS lyase HEAT-like repeat/i  ||
		$product =~ /^Sel1-like repeat$/i ||
      $product =~ /^Parallel beta-helix repeat$/i  ||
      $product =~ /^Helix-turn-helix motif$/i  ||
		$product =~ /^Transferase hexapeptide repeat$/i ||
      $product =~ /^Helix-turn-helix, type \d+$/i  ||
	   $product =~ /^(Thioredoxin|HTH|Potential Sulfotransferase) domain$/i ||
	   $product =~ /^S-layer domain$/i ||
      $product =~ /^(Thioredoxin|HTH) domain family$/i ||
      $product =~ /^Amino acid adenylation domain$/i ||
      $product =~ /^CopG-like DNA-binding$/i ||
      $product =~ /^Helix-turn-helix$/i ||
      $product =~ /^Helix-turn-helix HxlR type$/i ||
      $product =~ /^Alpha\/beta hydrolase fold$/i  ||
      $product =~ /^Hydrolase, alpha.beta fold/i  ||
      $product =~ /^Hydrolase, alpha.beta hydrolase fold/i  ||
      $product =~ /^unknown \w+-like protein$/i  ||
		$product =~ /^C-terminal half of Cry protein$/i ||
      $product =~ /^protein of unknown function/i ||
      $product =~ /^unknown protein [A-Z]\d+/i  ||
		$product =~ /^Similarities to phage protein/i ||
      $product =~ /^hypothetical membrane spanning protein$/i ||
      $product =~ /^Conserved hypothetical protein$/i ||
      $product =~ /^Uncharacterized conserved membrane protein/ ||
		$product =~ /^Similar to (hypothetical|bacteriophage) protein/i ||
		$product =~ /^Similar to short unknown protein/i ||
      $product =~ /^Truncated phage-like/ ||
		$product =~ /^Hypothetical phage protein$/i ||
      $product =~ /^Phage-related protein/i ||
		$product =~ /^Hypothetical transcriptional regulator$/i ||
		$product =~ /^Mitochondrial transporter/i ||
		$product =~ /^Uncharacterized low-complexity/ ||
      $product =~ /^Conserved hypothetical membrane protein$/i ||
      $product =~ /^integral membrane protein TIGR\d+$/i ||
	   $product =~ /^hypothetical two domain protein$/i ||
      $product =~ /^Conserved hypothetical integral membrane protein/ ||
      $product =~ /^[a-z]+ \(Conserved protein [a-z]+\)$/i ||
      $product =~ /^[a-z]+ \(Conserved membrane protein [a-z]+\)$/i ||
	   $product =~ /^Membrane protein, putative$/i ||
	   $product =~ /^Conserved protein\s*$/i ||
		$product =~ /Predicted transcriptional regulator with an addtional conserved domain/ ||
		$product =~ /transcriptional regulator with an additional conserved domain/i ||
		$product =~ /^Conserved (predicted|phage|domain|membrane|exported|hypothetical) protein/i ||
      $product =~ /^membrane spann?ing protein$/ ||
	   $product =~ /^inner membrane protein$/i ||
	   $product =~ /^Predicted small secreted protein$/i ||
	   $product =~ /^Uncharacterized conserved small protein/i ||
	   $product =~ /^Uncharacterized conserved protein (UCP|CAC)\d+/i ||
      $product =~ /^Outer membrane protein [A-Z]+\d+$/i ||
	   $product =~ /^Predicted protein\s*$/i ||
	   $product =~ /^Uncharacterized conserved$/i ||
	   $product =~ /putative uncharacterized protein/i ||
	   $product =~ /^Hypothetical (conserved|exported) protein$/i ||
	   $product =~ /^Hypothetical protein, partial$/i ||
	   $product =~ /^hypothetical protein membrane protein$/i ||
	   $product =~ /^Uncharacterized protein conserved in bacteria$/i ||
	   $product =~ /^Uncharacterized conserved$/i ||
	   $product =~ /^lipoprotein$/i ||
      $product =~ /^antigen$/i ||
	   $product =~ /^Bifunctional enzyme, contains/i ||
      $product =~ /^Group-specific protein$/i ||
	   $product =~ /^Extended ORF of/i ||
	   $product =~ /^Uncharacterized (membrane|conserved) protein$/i ||
	   $product =~ /^Exported membrane protein$/i ||
	   $product =~ /^Uncharacterized inner membrane protein$/i ||
	   $product =~ /^Uncharacterized conserved small protein-like protein/i ||
      $product =~ /^Predicted transcriptional regulator with an addtional conserved domain/i ||
      $product =~ /^Uncharacterized HTH-type transcriptional regulator/i ||
	   $product =~ /^Uncharacterized protein/i ||
      $product =~ /^uncharacterized domain \d+$/i ||
	   $product =~ /^Uncharacterized [.\d]+ kDa protein/i ||
	   $product =~ /^Possible (CDF|PET) family/ ||
	   $product =~ /^(Possible|Similarity|Uncharacterized|protein|toxin|predicted|hypothetical|exported)$/i || 
	   $product =~ /^\w+-related$/i ||
      $product =~ /^Enzyme\s?$/i ||
		$product =~ /^Transposase IS\d+.$/i ||
		$product =~ /^Bacillus cereus specific protein, uncharacterized$/i ||
		# product is 'hypothetical' if just an id with a vague name, e.g. 'Lin0391 protein'
      $product =~ /^Lin\d+ protein \(Lin\d+ protein\)$/ ||
		$product =~ /^[A-Z][a-z]{2}[A-Z]$/ ||
      $product =~ /^LP\d+G.\d+$/ ||
		$product =~ /^R.EcoHK31I protein$/i ||
      $product =~ /^lmo\d+$/i ||
		$product =~ /^Orf\d+\s*$/i ||
      $product =~ /^\d+orf\d+$/i ||
      $product =~ /^Orf\s+\d+[A-Z]+\s*$/i ||
      $product =~ /^Orf\s+[A-Z]+\d+\s*$/i ||
		$product =~ /^(IS66 Orf2|Orf2) like$/i ||
      $product =~ /^Orf\s+[A-Z]+\d+\s+putative/i ||
		$product =~ /^IS66 Orf2 like$/i ||
      $product =~ /^Orf\d+-like protein\s*$/i ||
      $product =~ /^Orf\d+-[A-Z]+\d+\s*$/i ||
      $product =~ /^Orf\S+\s+\(\S+\s+protein\)\s*$/i ||
      $product =~ /^Ig\s+hypothetical\s+\d+$/i ||
      $product =~ /^[a-z]{2}$/i ||
		$product =~ /^(HI|EF|Orf|Mob|MW|Blr|Cro|p|Blr|Orf-)\d+;?$/i ||
		$product =~ /^Possible \(\S+\) orf\d+/i ||
      $product =~ /^EF\d+\s+\(.+?\)$/i ||
      $product =~ /^Orf[A-Z]\s*$/ ||
		$product =~ /^GG\d+,?\s+isoform\s+[a-z]+$/i ||
      $product =~ /^\d+orf\d+\s*$/ ||
      $product =~ /^Cj\d+-\d+$/i ||
      $product =~ /^PXO\d+-\d+\s*$/i ||
      $product =~ /^UPI\S+\s+cluster$/i ||
      $product =~ /^Vng\d+[a-z]$/i ||
		$product =~ /^Zwa5b$/i ||
		$product =~ /^PclL$/i ||
		$product =~ /^\d+$/ ||
		$product =~ /^(Rep1|Doc|Vip2Ac|Zwa5A|Ugd|Sip1A|Vip1A\(BR\)|Tnp166|Blr5358|LAAC|MW2125)$/ ||
		$product =~ /^Gene \d+ protein$/ ||
		$product =~ /^Gp\d+\.\d+$/ ||
		$product =~ /^(19|Gne|Aec\d{1,2})$/ ||
		$product =~ /^protein\s+[\d{1,}a-z{1,}\/_-]+$/i ||
		$product =~ /^[\d{1,}a-z{1,}\/_-]+\s+protein$/i ||
      $product =~ /^lipoprotein\s+[\d{1,}a-z{1,}\/-_]+$/i ||
      $product =~ /^Maf-like protein\s+[\d{1,}a-z{1,}\/-_]+$/i ||
      $product =~ /^[\d{1,}a-z{1,}\/-_]+\s+protein\s+[\d{1,}a-z{1,}\/-_]+$/i 
					);

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
	return 3 if ( ($diff % 3) == 1 ); # confirmed
}

sub newFeatures {
	my ($self, @feat) = @_;

	push @{$self->{newFeatures}}, @feat if @feat;
	return @{$self->{newFeatures}};
}

sub run_tbl2asn {
	my ($self,$comment,$run) = @_;

	my ($tmplt, $outdir, $tbl2asn, $id) = 
	  ($self->template, $self->outdir, $self->executable, $self->id);

	if ( $run > 1) {
		for my $suffix ( qw( val sqn gbf ) ) {
			system "mv $outdir/$id.$suffix $outdir/$id.$suffix.orig" if ( -e "$outdir/$id.$suffix" );
		}
		system "mv discrp discrp.orig" if ( -e "discrp" );
	}

	if ( $outdir && $tmplt && $tbl2asn ) {
		my $cmd = "$tbl2asn -s -j [gcode=11] -V vb -Z discrp -t $tmplt.sbt -p $outdir -y \"$comment\"";
		print "tbl2asn command is \'$cmd\'\n" if $self->debug;
		`$cmd`;
		return 1;
	}
	die "Problem running tbl2asn";
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
sub get_dup_rnas {
	my $self = shift;
	my @todelete = ();
	my $id = $self->id;
	my $gene1;

	my @overlapgenes = $self->get_from_discrp('OVERLAPPING_GENES');

	$gene1 = shift @overlapgenes;

	while ( my $gene2 = shift @overlapgenes ) {

		if ( $gene1 =~ /_r\d+/ && $gene2 =~ /_r\d+/ ) {

			my ($g1start,$g1end) = $gene1 =~ /contig[\d.]+:c?(\d+)[<>]?-[<>]?(\d+)/;
			my ($g2start,$g2end) = $gene2 =~ /contig[\d.]+:c?(\d+)[<>]?-[<>]?(\d+)/;

			my ($g1ctg) = $gene1 =~ /(contig[.\d]+)/;
			my ($g2ctg) = $gene2 =~ /(contig[\d.]+)/;

			# remove identical rRNAs and rRNAs inside rRNAs
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

			# remove the smaller of 2 overlapping rRNAs
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

	while ( my $gene2 = shift @overlapgenes ) {

		# Note that this loop does not collect rRNAs or tRNAs

		my ($g1start,$g1end) = $gene1 =~ /contig[\d.]+:c?(\d+)[<>]?-[<>]?(\d+)/;
		my ($g2start,$g2end) = $gene2 =~ /contig[\d.]+:c?(\d+)[<>]?-[<>]?(\d+)/;

		my $g1len = abs($g1start - $g1end);
		my $g2len = abs($g2start - $g2end);

		my ($g1ctg) = $gene1 =~ /(contig[.\d]+)/;
		my ($g2ctg) = $gene2 =~ /(contig[\d.]+)/;

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
sub get_rna_overlaps {
	my $self = shift;
	my @todelete = ();
	my $gene1;
	my $id = $self->id;

	my @overlaps = $self->get_from_discrp('RNA_CDS_OVERLAP');

	$gene1 = shift @overlaps;

	while ( my $gene2 = shift @overlaps ) {

      my ($g1start,$g1end) = $gene1 =~ /contig[\d.]+:c?(\d+)[<>]?-[<>]?(\d+)/;
      my ($g2start,$g2end) = $gene2 =~ /contig[\d.]+:c?(\d+)[<>]?-[<>]?(\d+)/;
      my ($g1ctg) = $gene1 =~ /(contig[\d.]+)/;
      my ($g2ctg) = $gene2 =~ /(contig[\d.]+)/;

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
		$readflag = 1 if ( /^DiscRep_SUB:$header/ );
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

	system "mv $dir/$id.tbl $dir/$id.tbl.orig" if ( -e "$dir/$id.tbl" );

	open MYOUT,">$dir/$id.tbl" or die "Cannot write to $dir/$id.tbl";
	print "Writing to $dir/$id.tbl\n" if $self->debug;

	for my $contig ( @{$tbl} ) {

		print MYOUT ">Features " . $contig->{contigname} . "\n";

	 FEAT:
		for my $feat ( sort sort_by_loc keys %{$contig} ) {
			next FEAT if ( $feat eq 'contigname' );
			print MYOUT $feat;
			print MYOUT $contig->{$feat};
		}
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

	# transform array of lines into hash of names plus locations
	# Gene    yberc_r90       lcl|contig00186:128269->128665  yberc_r90
	# CDS     hypothetical protein    lcl|contig00890:c203-<1 yberc_39200
	for my $delete (@todelete) {
		$delete =~ /^\S+\s+[^|]+\|(contig[\d.]+):c?(\d+)[<>]?-[<>]?(\d+)/;
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

	#unlink "$template.sbt" if ( -e "$template.sbt" );

	unlink "discrp.orig" if -e "discrp.orig";
	system "mv discrp $dir" if -e "discrp";

	system "mv $id.agp $dir" if -e "$id.agp";

	for my $suffix ( qw(gbf val tbl sqn) ) {
		unlink "$dir/$id.$suffix.orig" if -e "$dir/$id.$suffix.orig";
	}

}

sub outdir {
	my ($self,$id) = @_;

	if ($id) {
		# make output directory
		my $outdir = "$id-gbsubmit.out.d";
		`rm -rf $outdir; mkdir -p $outdir`;
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
	my $strain = $self->strain;

	local $/ = undef;

   my $template = $self->template;
   my $asn = $template . ".sbt";
   unlink $asn if (-e $asn);

	open MYIN,$template or die "Cannot open file $template";
	my $text = <MYIN>;
	$text =~ s/<BACTERIALSPECIES>/$species/;
   $text =~ s/<BACTERIALSTRAIN>/$strain/;

	open MYOUT,">$asn" or die "Cannot write to file $asn";
	print MYOUT $text;

}

sub make_top_comment {
	my $self = shift;
	my ($totalLen,$totalReads) = 0;

	my $readsPerBase = $self->readsPerBase;
	my $comment = '';

	if ( $readsPerBase ) {
		for my $len ( keys %{$readsPerBase} ) {
			$totalLen += $len;
			$totalReads += ( $len * $readsPerBase->{$len} );
		}

		$comment = "Genome coverage: " . 
		  sprintf("%.1f",($totalReads/$totalLen)) . "X";
	}

	if ( $self->strain !~ /ATCC/ ) {
		if ( $comment ) {
			$comment = "Bacteria available from BDRD~" . $comment;
		} else {
			$comment = "Bacteria available from BDRD";
		}
	}

	$comment;
}

sub create_qual {
	my ($self,$qualfile) = @_;

	my $contignames = $self->contigs;
	my $outdir = $self->outdir;
	my $id = $self->id;
	my %qualities = ();

	if ( -e $qualfile ) {
		my $in = Bio::SeqIO->new(-file => $qualfile,
										 -format => "qual" );

		while ( my $qual = $in->next_seq ) {
			my $primary_id = $qual->primary_id;
			$qualities{$primary_id} = $qual;
		}

		my $out = Bio::SeqIO->new(-file => ">>$outdir/$id.qvl",
										  -format => "qual" );

		# A name could be something like "contig01112.6" or "contig01112"
		for my $name ( @{$contignames} ) {

			for my $key ( keys %qualities ) {
				if ( $name =~ /$key/ ) {
					my $obj = $qualities{$key};
					$obj->display_id($name);
					$out->write_seq($obj)
				}
			}
		}

	} else {
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
	my ($self,$gbk) = @_;
	my $taxid = $self->taxid or die "No taxid found";
	my $date = $self->get_date;
	my $id = $self->id;
   my $organism = $self->organism;
   my $strain = $self->strain;
	my %frames = ('1','+','-1','-');

	# hack!
	#my @gbks = <*rnammer.out.gbk>;
	#die "No *rnammer.out.gbk file found, can not create an *.agp file" unless ( -e $gbks[0] );

	open MYAPG,">>$id.agp" or die "Cannot create file $id.agp";

	my $in = Bio::SeqIO->new( -file => $gbk, -format => 'genbank' );
	my $seq = $in->next_seq;
	my @contigs = grep { $_->primary_tag eq 'fasta_record' } $seq->get_SeqFeatures;

my $text = "# ORGANISM: $organism
# TAX_ID: $taxid
# ASSEMBLY NAME: $id
# ASSEMBLY DATE: $date
# GENOME CENTER: NMRC
# DESCRIPTION: $organism $strain chromosome, whole genome shotgun
";

  print MYAPG $text;

   my $count = 1;
   my $pos = 1;
   my $spacer_len = 100;

	for my $contig (@contigs) {
	  my @names = $contig->get_tag_values('name');
	  my $len = ($contig->end) - ($contig->start);
	  print MYAPG "Chr1\t" . $pos . "\t" . ($len + $pos) . "\t" . $count . "\tW\t" . $names[0] . "\t1\t" . 
		 ($len + 1) . "\t" . ($frames{$contig->strand}) . "\n";
	  $pos = $len + $pos + 1;
	  $count++;
	  print MYAPG "Chr1\t" . $pos .  "\t" . ($pos + $spacer_len - 1) . "\t" . $count . "\tN\t$spacer_len\tfragment\tyes\n";
	  $pos += $spacer_len;
	  $count++;
  }


}

sub get_date {
	my $self = shift;

	use Date::Format;
	return time2str("%d-%B-%Y",time);
}

sub edit_definition {
	my ($self,$desc) = @_;

	# Edit DEFINITION line, input looks something like: 
	# [organism=Yersinia bercovieri] [strain=ATCC_43970] [gcode=11] [date=7-3-2008]
	my ($organism) = $desc =~ /\[organism=([^]]+)/;
	my ($strain) = $desc =~ /\[strain=([^]]+)/;
	my ($gcode) = $desc =~ /\[gcode=([^]]+)/;
	my $definition = '';

	if ($strain eq undef) {
		$strain = '';
	} else {
		$strain =~ s/_/ /g;
	}

	$self->organism($organism);
	$self->strain($strain);

	if ( $strain =~ /ATCC/ ) {
		my ($id) = $strain =~ /ATCC\s+(.+)/;
      $definition = "[organism=$organism][strain=$strain][culture_collection=ATCC:$id][gcode=$gcode][tech=wgs]";
	} else {
		$definition = "[organism=$organism][strain=$strain][gcode=$gcode][tech=wgs]";
	}
	$definition;
}

# get names, see if any are duplicated - /name="contig01098"
#     fasta_record    1..529
#                     /name="contig00231"
#                     /note="untiled"
#                     /note="Calculated sequence coverage = 5.36 reads per base"
sub make_namemap {
	my ($self,$file) = @_;
	my $namemap;

	my $in = Bio::SeqIO->new(-file => $file, -format => 'genbank');
	my $seq = $in->next_seq;
	for my $feat ($seq->get_SeqFeatures) {
		if ( $feat->primary_tag eq 'fasta_record' ) {
			my @names = $feat->get_tag_values('name');
			$namemap->{$names[0]}->{num}++;
	      # $namemap{$1}++ if ( /\/name="([^"]+)/ );
			$namemap->{$names[0]}->{len} = ($feat->end) - ($feat->start);
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


sub lastBase {
	my ($self,$base) = @_;

	if ($base) {
		$self->{lastBase} = $base;
		return $self->{lastBase};
	} else {
		return $self->{lastBase};
	}
}

sub executable {
	my ($self,$exec) = @_;

	if ($exec) {
		$self->{executable} = $exec;
		return $self->{executable};
	} else {
		return $self->{executable};
	}
}

sub namemap {
   my ($self,$ref) = @_;

   if ($ref) {
      $self->{namemap} = $ref;
      return $self->{namemap};
   } else {
      return $self->{namemap};
   }
}

sub cutoff {
   my ($self,$cut) = @_;

   if ($cut) {
      $self->{cutoff} = $cut;
      return $self->{cutoff};
   } else {
		return 0 if (! defined $self->{cutoff} );
      return $self->{cutoff};
   }
}

sub debug {
	my ($self,$debug) = @_;

	if ($debug) {
		$self->{debug} = $debug;
		return $self->{debug};
	} else {
		return $self->{debug};
	}
}

sub template {
	my ($self,$template) = @_;

	if ($template) {
		$self->{template} = $template;
		return $self->{template};
	} else {
		return $self->{template};
	}
}

sub taxid {
   my ($self,$id) = @_;

   if ($id) {
      $self->{taxid} = $id;
      return $self->{taxid};
   } else {
      return $self->{taxid};
   }
}

sub organism {
   my ($self,$id) = @_;

   if ($id) {
      $self->{organism} = $id;
      return $self->{organism};
   } else {
      return $self->{organism};
   }
}

sub id {
	my ($self,$id) = @_;

	if ($id) {
		$self->{id} = $id;
		return $self->{id};
	} else {
		return $self->{id};
	}
}

sub strain {
	my ($self,$strain) = @_;
	$self->{strain} = $strain if ($strain);
	$self->{strain};
}

sub contigs {
	my ($self,$contigs) = @_;
	$self->{contigs} = $contigs if ($contigs);
	$self->{contigs};
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


