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

MARC::CRT

=head1 SYNOPSIS

Do not use this module directly. This module is used by diya.pm when it
runs a pipeline.

=head1 DESCRIPTION

Template Perl module. A parser() method is required.

=head1 AUTHOR

Brian Osborne, briano@bioteam.net

=cut

# add the new module name here
package diya::MARC::CRT;

use strict;
# simplest approach
use base 'diya';

# alternate approach, if a variable needs to be imported
# use vars qw(@ISA);
# use diya qw();
# @ISA = qw(diya);


=head2 parse

 Name    : parse
 Usage   : $parser->parse($diya)
 Function: parse a program output file - every parser module must have a 
           parse() method
 Returns : 1 on success
 Args    : a diya object
 Example : 

=cut

sub parse {
	my ($self,$diya) = @_;

	my $LOCUS_TAG_NUMBER = 0;

	# Parse CRT output, get features back
	my $out = $diya->_outputfile('MARC::CRT');
	print "Parsing " . $out . "\n" if $diya->verbose;
	my @crisprs = parse_crt($out);
	print "Found CRISPRS\n" if ( $diya->verbose && @crisprs );

	my $gbkin = $diya->_outputfile("MARC::rnammer");
	my $seqin = Bio::SeqIO->new(-file => "$gbkin.gbk", -format => 'genbank');
	my $seq = $seqin->next_seq;

	# Add any new features
	for my $crispr ( @crisprs ) {
		my %tag;
		$tag{locus_tag} = $seq->display_id . "_c" . ($LOCUS_TAG_NUMBER += 10);
		$crispr->set_attributes(-tag => \%tag);
		$crispr->source_tag('CRISPR Recognition Tool');
		$seq->add_SeqFeature($crispr);
	}

	# Sort features by location
	my @features = $seq->remove_SeqFeatures;
	@features = sort { $a->start <=> $b->start } @features;
	$seq->add_SeqFeature(@features);

	# Output to Genbank file
	my $gbkout = $out . ".gbk";
	my $seqout = Bio::SeqIO->new(-format => 'genbank',
				     -file   => ">$gbkout");
	$seqout->write_seq($seq);

}

sub parse_crt {
	my $file = shift;
	my $txt = read_file($file);
	my @features;

	# CRISPR 1   Range: 123524 - 124955
	while ( $txt =~ /CRISPR\s+\d+\s+Range:\s+(\d+)\s+-\s+(\d+)/g ) {

		my $feat = new Bio::SeqFeature::Generic(-start       => $1,
    		                                    -end         => $2,
        		                                -strand      => 1,
							-note => 'CRISPR',
            		                            -primary_tag => 'repeat_region');
        push @features,$feat;
    }

	return 0 if ( ! @features );

	@features;
}

sub read_file {
	my $file = shift;
	local $/ = undef;
	open MYIN,$file;
	my $txt = <MYIN>;
	$txt;
}

1;

__END__

