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

package diya::BDRD::rnammer;

use strict;
use base 'diya';
use Bio::Tools::GFF;
use Bio::SeqIO;
use Bio::SeqFeature::Generic;

sub parse {
    my ( $self, $diya ) = @_;

    my $LOCUS_TAG_NUMBER = 0;

    #my $gbk = $diya->_sequence;

    my $rin = $diya->_outputfile('BDRD::tRNAscanSE') . ".gbk";
    my $in  = Bio::SeqIO->new( -file => $rin, -format => 'genbank' );
    my $gbk = $in->next_seq;

    my $rout = $diya->_outputfile('BDRD::rnammer');

    open MYIN, $rout or die "Cannot read file $rout";

    while ( my $line = <MYIN> ) {

        if ( $line =~
/^\S+\s+(\w+)-([\d.]+)\s+rRNA\s+(\d+)\s+(\d+)\s+(\S+)\s+[+-]\s+\.\s+(\d+s)/
          )
        {
            my %tags;

            $tags{locus_tag} =
              $gbk->display_id . "_r" . ( $LOCUS_TAG_NUMBER += 10 );
            $tags{note}      = "$1 score=$5";
            $tags{product}   = "$6 ribosomal rna";
            $tags{inference} = "profile:$1:$2";

            my $feat = Bio::SeqFeature::Generic->new(
                -primary    => 'rRNA',
                -source_tag => 'rnammer',
                -start      => $3,
                -end        => $4,
                -tag        => {%tags}
            );
            $gbk->add_SeqFeature($feat);
        }
    }

    # Sort features by location
    my @features = $gbk->remove_SeqFeatures;

    # sort features by start position
    @features = sort { $a->start <=> $b->start } @features;
    $gbk->add_SeqFeature(@features);

    my $seqo = Bio::SeqIO->new(
        -format => 'genbank',
        -file   => ">$rout.gbk"
    );
    $seqo->write_seq($gbk);
}

1;

__END__


##gff-version2
##source-version RNAmmer-1.2
##date 2008-11-11
##Type DNA
# seqname           source                      feature     start      end   score   +/-  frame  attribute
# ---------------------------------------------------------------------------------------------------------
NC_004061-contig5       RNAmmer-1.2     rRNA    482     2021    1855.8  +       .       16s_rRNA        
# ---------------------------------------------------------------------------------------------------------

    gene           130..219
                    /locus_tag="yberc_r10"
    rRNA           130..219
                    /locus_tag="yberc_r10"
                    /product="ribosomal RNA"
                    /inference="profile:RNAMMER:1.2"
                    /note="RNAMMER score=62.9"
