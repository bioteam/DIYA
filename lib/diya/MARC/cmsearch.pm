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

cmsearch

=head1 SYNOPSIS

Do not use this module directly. This module is used by diya.pm when it
runs a pipeline.

=head1 DESCRIPTION

Template Perl module. A parser() method is required.

=head1 AUTHOR

Brian Osborne, briano@bioteam.net

=cut

# add the new module name here
package diya::MARC::cmsearch;

use strict;
use Bio::SeqIO;
use Bio::SeqFeature::Generic;
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
	my @rnas;

	# Parse cmsearch output, get features back
	my $outputfile = $diya->_outputfile('MARC::cmsearch');
	print "Parsing $outputfile\n" if $diya->verbose;
	@rnas = parse_cmsearch($outputfile);
	print "Found " . scalar @rnas . " sRNAs\n" if ( $diya->verbose && @rnas );

	# Get the annotated sequence from the previous step
	my $gbkin = $diya->_outputfile("MARC::phobos");
	my $seqin = Bio::SeqIO->new(-file => "$gbkin.gbk", -format => 'genbank');
	my $seq = $seqin->next_seq;

	# Add any new features from cmsearch
	for my $rna ( @rnas ) {
		# $rna->source_tag('cmsearch');
		$seq->add_SeqFeature($rna);
	}

	# Sort features by location
	my @features = $seq->remove_SeqFeatures;
	@features = sort { $a->start <=> $b->start } @features;
	$seq->add_SeqFeature(@features);

	# Output to Genbank file
	my $gbkout = $outputfile . ".gbk";
	my $seqout = Bio::SeqIO->new(-format => 'genbank',
				                 -file   => ">$gbkout");
	$seqout->write_seq($seq);

}

sub parse_cmsearch {
	my $file = shift;
	my @features;

# Example:
# Spot_42     gi|315134697|dbj|AP012030.1|     4028994     4029112      1    118    124.06  5.54e-29   43
# Spot_42     gi|315134697|dbj|AP012030.1|      284360      284447      1    118     26.13  1.85e-03   40

	open MYIN,$file;
	while ( <MYIN> ) {

		/^\s+(\S+)\s+\S+\s+(\d+)\s+(\d+)\s+\d+\s+\d+\s+[\d\.]+\s+/;

		if ( $1 && $2 && $3 ) {

			my $feat = new Bio::SeqFeature::Generic(
					     -start  => $2,
    		             -end    => $3,
        		         -strand => 1,
        		         -tag    => { note        => "$1 sRNA",
        		                      ncRNA_class => 'other',
        		                      inference   => "profile:cmsearch:1.0.1" },
            		     -primary_tag => 'ncRNA');

        	push @features,$feat;
        	print "Found $1 sRNA at $2 and $3\n";
    	}
    }

   @features;
}

1;

__END__
