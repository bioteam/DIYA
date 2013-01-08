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

tRNAscanSE.pm

=head1 SYNOPSIS

Do not use this module directly. This module is used by diya.pm when it
runs a pipeline.

=head1 DESCRIPTION

A tRNA finding component of the diya project which serves
as a wrapper and result parser for tRNAscan-SE.

=head1 AUTHOR

Andrew Stewart, andrew.stewart@med.navy.mil

=head1 CONTRIBUTORS

Tim Read, timothy.read@med.navy.mil
Brian Osborne, briano@bioteam.net

=cut

package diya::BDRD::tRNAscanSE;

use strict;
use base 'diya';
use Bio::Tools::tRNAscanSE;

sub parse {
    my ( $self, $diya ) = @_;

    my $LOCUS_TAG_NUMBER = 0;

    my $out = $diya->_outputfile('BDRD::tRNAscanSE');
    print "Parsing " . $out . "\n" if $diya->verbose;

    # Parse tRNAscan output
    my $parser = Bio::Tools::tRNAscanSE->new(
        -file    => "$out",
        -genetag => 'tRNA'
    );

    my $gbk = $diya->_outputfile("BDRD::rpsblastCDS");
    my $in  = Bio::SeqIO->new( -file => "$gbk.gbk", -format => 'genbank' );
    my $seq = $in->next_seq;

    # Parse the results
    while ( my $feat = $parser->next_prediction ) {
        my %tag;
        $tag{locus_tag} = $seq->display_id . "_t" . ( $LOCUS_TAG_NUMBER += 10 );
        $feat->set_attributes( -tag => \%tag );
        $feat->add_tag_value( 'inference', 'profile:tRNAscan-SE:1.23' );
        $seq->add_SeqFeature($feat);
    }

    # Sort features by location
    my @features = $seq->remove_SeqFeatures;
    @features = sort { $a->start <=> $b->start } @features;
    $seq->add_SeqFeature(@features);

    # Output
    my $outfile = $out . ".gbk";
    my $seqo    = Bio::SeqIO->new(
        -format => 'genbank',
        -file   => ">$outfile"
    );
    $seqo->write_seq($seq);
}

1;

__END__




