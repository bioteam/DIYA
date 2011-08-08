# $Id: tRNAscanSE.pm 297 2008-12-18 15:36:01Z briano $
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

package diya::PipeLineOne::phagefinder;

use strict;
use base 'diya'; 

sub parse {
	my ($self,$diya) = @_;

	my $LOCUS_TAG_NUMBER = 0;

	my $out = $diya->_outputfile('PipeLineOne::phagefinder');
	print "Parsing " . $out . "\n" if $diya->verbose;

	# Retrieve Phage_Finder output
#	my $parser = Bio::Tools::tRNAscanSE->new(-file => "$out",
#								             -genetag	=> 'tRNA');

	open INFILE, $out;
	my @tabfile = <INFILE>;
#	my $outdir = $MYSEQID . '_dir';
#	my @tabfile = `cat $outdir/*.tab`;

	my $gbk = $diya->_outputfile("PipeLineOne::rnammer"); # TODO
	my $in = Bio::SeqIO->new(
		-file => "$gbk.gbk", 
		-format => 'genbank',
	);
	my $seq = $in->next_seq;

	# my $seq = $diya->_sequence;

	# parse the results
	for my $line (@tabfile) {
	
		my @col = split "\t", $line;
	
		my ($start, $end, $strand);
		if ($col[3] <= $col[4]) {
			$start = $col[3];
			$end = $col[4]; 
			$strand = '+';
		} else {
			$start = $col[4];
			$end = $col[3];
			$strand = '-';
		}
	
		my %tag;		
		$tag{locus_tag} = $seq->display_id . "_p" . ($LOCUS_TAG_NUMBER += 10);
		my $feat = Bio::SeqFeature::Generic->new( 
            -start        => $start, 
            -end          => $end,
            -strand       => $strand,
            -primary      => $col[7], # -primary_tag is a synonym
            -source_tag   => 'Phage_Finder',
            -display_name => "$col[6] $col[7]",
        );
		$feat->set_attributes(-tag => \%tag);
		$seq->add_SeqFeature($feat);
	}

	# Sort features by location
	my @features = $seq->remove_SeqFeatures;
	# sort features by start position
	@features = sort { $a->start <=> $b->start } @features;
	$seq->add_SeqFeature(@features);

	# Output

	my $outfile = $out . ".gbk";

	my $seqo = Bio::SeqIO->new(-format	=> 'genbank',
							   -file	=> ">$outfile");
	$seqo->write_seq($seq);

}

1;

__END__



my $INFILE;
my $OUTFILE;
my $CONF;
my $setup;
my $FLAG_B = 1;
my $CLEANUP;
my $VERBOSE;
my $LOCUS_TAG_NUMBER = 0;
my $BIN;

GetOptions(	"conf=s"		=> \$CONF,
			"outfile=s"		=> \$OUTFILE,
			"bin=s"			=> \$BIN,
			"cleanup|c"		=> \$CLEANUP,
			"verbose|v"		=> \$VERBOSE,
			"bacterial"		=> \$FLAG_B,
);

# import configuration directly
if ($CONF && -e $CONF) {
	$setup = XMLin($CONF);
}

$INFILE = shift @ARGV;
$OUTFILE = $INFILE unless ($OUTFILE);
my $fileroot = $INFILE;
$fileroot =~ s/\.gbk//;

if($BIN) {
	$BIN .= '/' unless ($BIN =~ /$\//);
}

#-------------------------------------------------------------------------------
# INPUT
# Input genbank file
my $seqi = Bio::SeqIO->new(
	-format	=> 'genbank',
	-file	=> $INFILE);
my $seq = $seqi->next_seq;

# Convert to a temporary fasta file
unless (-e "$fileroot.fsa") {
	my $temp_fasta = Bio::SeqIO->new(
		-format => 'fasta',
		-file	=> ">$fileroot.fsa");
	$temp_fasta->write_seq($seq);
}

#-------------------------------------------------------------------------------
# Execute tRNAscan on temp fasta file with parameters
my $command = $BIN || $setup->{conf}->{trnascan}->{bin} || "tRNAscan-SE";
my @params = (
	'-B',
	-o	=> "$fileroot.trnascan",
	" $fileroot.fsa",
);


