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

sRNAscanner

=head1 SYNOPSIS

Do not use this module directly. This module is used by diya.pm when it
runs a pipeline.

=head1 DESCRIPTION

Template Perl module. A parser() method is required.

=head1 AUTHOR

Brian Osborne, briano@bioteam.net

=cut

# add the new module name here
package diya::MARC::sRNAscanner;

use strict;
use Bio::SeqIO;
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

	my $srnaoutput = '/usr/local/share/apps/sRNAscanner/output_files/sRNA.txt';

	# Parse sRNAscanner output, get features back
	print "Parsing $srnaoutput\n" if $diya->verbose;
	my @rnas = parse_sRNAscanner($srnaoutput);
	print "Found " . scalar @rnas . " sRNAs\n" if ( $diya->verbose && @rnas );

	# Get the annotated sequence from the previous step
	my $gbkin = $diya->_outputfile("MARC::phobos");
	my $seqin = Bio::SeqIO->new(-file => "$gbkin.gbk", -format => 'genbank');
	my $seq = $seqin->next_seq;

	# Add any new features from sRNAscanner
	for my $rna ( @rnas ) {
		my %tag;
		$tag{note} = 'sRNA';
		$tag{ncRNA_class} = 'other';
		$rna->set_attributes(-tag => \%tag);
		$rna->source_tag('sRNAscanner');
		$seq->add_SeqFeature($rna);
	}

	# Sort features by location
	my @features = $seq->remove_SeqFeatures;
	@features = sort { $a->start <=> $b->start } @features;
	$seq->add_SeqFeature(@features);

	# Output to Genbank file
	my $out = $diya->_outputfile('MARC::sRNAscanner');
	my $gbkout = $out . ".gbk";
	my $seqout = Bio::SeqIO->new(-format => 'genbank',
				     -file   => ">$gbkout");
	$seqout->write_seq($seq);

}

sub parse_sRNAscanner {
	my $file = shift;
	my @features;

# >sRNA|39708|39888|
# TTATGTCTCTGTTGTAAAAGTCACACCGGATAGCATGAAATTAATGAAACTTCGAATGGGAATAATCTC

	my $in = Bio::SeqIO->new(-file => $file, -format => 'fasta');

	while ( my $seq = $in->next_seq ) {

		$seq->desc =~ /sRNA\|(\d+)\|(\d+)/;

		my $feat = new Bio::SeqFeature::Generic(-start       => $1,
    		                                    -end         => $2,
        		                                -strand      => 1,
            		                            -primary_tag => 'ncRNA');
        push @features,$feat;
        print "Found sRNA at $1 and $2\n";
    }

	return 0 if ( ! @features );

	@features;
}

1;

__END__
