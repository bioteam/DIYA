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

phobos

=head1 SYNOPSIS

Do not use this module directly. This module is used by diya.pm when it
runs a pipeline.

=head1 DESCRIPTION

Template Perl module. A parser() method is required.

=head1 AUTHOR

Brian Osborne, briano@bioteam.net

=cut

# add the new module name here
package diya::MARC::phobos;

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

	# Parse phobos output, get features back
	my $out = $diya->_outputfile('MARC::phobos');
	print "Parsing " . $out . "\n" if $diya->verbose;
	my @repeats = parse_phobos($out);
	print "Found repeats\n" if ( $diya->verbose && @repeats );

	# Get the annotated sequence from the previous step
	my $gbkin = $diya->_outputfile("MARC::CRT");
	my $seqin = Bio::SeqIO->new(-file => "$gbkin.gbk", -format => 'genbank');
	my $seq = $seqin->next_seq;

	# Add any new features from phobos
	for my $repeat ( @repeats ) {
		my %tag;
		$tag{rpt_type} = 'tandem';
		$repeat->set_attributes(-tag => \%tag);
		$repeat->source_tag('phobos');
		$seq->add_SeqFeature($repeat);
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

sub parse_phobos {
	my $file = shift;
	my $txt = read_file($file);
	my @features;

# contig00007	Phobos	tandem-repeat	8849	8866	100.00	.	.	Name="repeat_region 8849-8
# 866 unit_size 9 repeat_number 2.000 perfection 100.000 unit ATCGCCGCC"

	while ( $txt =~ /^\S+\s+Phobos\s+tandem-repeat\s+(\d+)\s+(\d+)/g ) {

		my $feat = new Bio::SeqFeature::Generic(-start       => $1,
    		                                    -end         => $2,
        		                                -strand      => 1,
        		                                -tag         => {inference => 'phobos' },
            		                            -primary_tag => 'repeat_region');
        push @features,$feat;
    }

	return 0 if ( ! @features );

	@features;
}

sub read_file {
	my $file = shift;
	my $txt;

	open MYIN,$file;
	while (<MYIN>) {
		$txt .= $_ if ( ! /^#/ );	
	}

	$txt;
}

1;

__END__
