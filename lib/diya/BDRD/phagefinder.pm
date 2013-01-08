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

phagefinder.pm

=head1 SYNOPSIS

Do not use this module directly. This module is used by diya.pm when it
runs a pipeline.

=head1 DESCRIPTION

A phage finding component of the diya project which serves
as a wrapper and result parser for PhageFinder.

=head1 AUTHOR

Mohit Patel

=head1 CONTRIBUTORS

Brian Osborne, briano@bioteam.net

=cut

package diya::BDRD::phagefinder;

use strict;
use base 'diya';

sub parse {
    my ( $self, $diya ) = @_;

    my $LOCUS_TAG_NUMBER = 0;

    my $out = $diya->_outputfile('BDRD::phagefinder');
    print "Parsing " . $out . "\n" if $diya->verbose;

    my $gbk = $diya->_outputfile("BDRD::cmsearch");
    my $in  = Bio::SeqIO->new(
        -file   => "$gbk.gbk",
        -format => 'genbank',
    );
    my $seq = $in->next_seq;

    # Retrieve Phage_Finder output
    open INFILE, $out;
    my @tabfile = <INFILE>;

    for my $line (@tabfile) {

        my @col = split "\t", $line;

        my ( $start, $end, $strand );
        if ( $col[3] <= $col[4] ) {
            $start  = $col[3];
            $end    = $col[4];
            $strand = '+';
        }
        else {
            $start  = $col[4];
            $end    = $col[3];
            $strand = '-';
        }

        my %tag;
        $tag{locus_tag} = $seq->display_id . "_p" . ( $LOCUS_TAG_NUMBER += 10 );
        my $feat = Bio::SeqFeature::Generic->new(
            -start        => $start,
            -end          => $end,
            -strand       => $strand,
            -tag          => { inference => "profile:Phage_Finder:1.0" },
            -primary_tag  => 'misc_feature',
            -note         => 'bacteriophage',
            -source_tag   => 'Phage_Finder',
            -display_name => $col[6] . ' ' . $col[7],
        );
        $feat->set_attributes( -tag => \%tag );
        $seq->add_SeqFeature($feat);
    }

    # Sort features by location
    my @features = $seq->remove_SeqFeatures;

    # sort features by start position
    @features = sort { $a->start <=> $b->start } @features;
    $seq->add_SeqFeature(@features);

    # Output
    my $outfile = $out . ".gbk";

    my $seqo = Bio::SeqIO->new(
        -format => 'genbank',
        -file   => ">$outfile"
    );
    $seqo->write_seq($seq);

}

1;

__END__
