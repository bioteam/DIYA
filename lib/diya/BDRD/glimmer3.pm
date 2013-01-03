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

glimmer3.pm

=head1 SYNOPSIS

Do not use this module directly. This module is used by diya.pm when it
runs a pipeline.

=head1 DESCRIPTION

The gene prediction component of the Do It Yourself Annotator utilizing glimmer
for ab initio detection of potential coding sequences.  

This script has the glimmer3 package in mind, so use with older versions 
of glimmer may require some tampering.  Speficially, use of shell scripts
in glimmer3 to control the individual processes (longorfs, etc) has been 
used here, as opposed to running the individual glimmer programs.

=head1 AUTHOR

Andrew Stewart, andrew.stewart@med.navy.mil

=head1 CONTRIBUTORS

Tim Read, timothy.read@med.navy.mil
Brian Osborne, briano@bioteam.net

=cut

package diya::BDRD::glimmer3;

use strict;
use vars qw(@ISA);
use diya qw($MYSEQID);
@ISA = qw(diya);

use Bio::Tools::Glimmer;
use Bio::SeqIO;
use Data::Dumper;

sub parse {
    my ( $self, $diya ) = @_;

    my $PRIMARY_TAG      = "gene";
    my $LOCUS_TAG_NUMBER = 0;        # starting number
    my ( $seq, $out );

    if ( defined $MYSEQID ) {
        $out = $MYSEQID;
    }
    else {
        $out = $diya->_outputfile('BDRD::glimmer3');
    }
    print "Parsing $out.predict\n" if $diya->verbose;

    my $parser = Bio::Tools::Glimmer->new(
        -file   => "$out.predict",
        -format => 'Glimmer'
    );

    if ( defined $MYSEQID ) {
        my $gbk = "$MYSEQID.gbk";

        if ( !-e $gbk ) {
            my $dir = $diya->outputdir;
            $gbk = "$dir/$gbk";
        }

        my $in = Bio::SeqIO->new(
            -file   => "$gbk",
            -format => 'genbank'
        );
        $seq = $in->next_seq;

        print "Sequence is $gbk\n" if $diya->verbose;

    }
    else {
        $seq = $diya->_sequence;
    }

    while ( my $feature = $parser->next_feature ) {
        my %tags;
        $tags{locus_tag} = $MYSEQID . '_' . ( $LOCUS_TAG_NUMBER += 10 );
        $tags{inference} = 'ab initio prediction:glimmer3:3.0.2';

        my $feat = Bio::SeqFeature::Generic->new(
            -primary    => $PRIMARY_TAG,
            -source_tag => 'glimmer3',
            -start      => $feature->start,
            -end        => $feature->end,
            -strand     => $feature->strand,
            -tag        => {%tags}
        );

        # Add feature to seq object
        $seq->add_SeqFeature($feat);
    }

    # Sort features by location
    my @features = $seq->get_SeqFeatures;

    # sort features by start position
    @features = sort { $a->start <=> $b->start } @features;
    $seq->remove_SeqFeatures;
    $seq->add_SeqFeature(@features);

    # Output
    my $outfile = $diya->_outputfile('BDRD::glimmer3');

    my $seqo = Bio::SeqIO->new(
        -format => 'genbank',
        -file   => ">$outfile.gbk"
    );
    $seqo->write_seq($seq);
    $diya->_sequence($seq);

}

1;

__END__



my $setup;
$GLIMMER = "glimmer3";
$BIN;

GetOptions(	"conf=s"				=> \$CONF,
			"outfile=s"				=> \$OUTFILE,
			"locus_tag_prefix=s"	=> \$LOCUS_TAG_PREFIX,
			"source_tag=s"			=> \$SOURCE_TAG,
			"primary_tag=s"			=> \$PRIMARY_TAG,
			"glimmer=s"				=> \$GLIMMER,
			"verbose"				=> \$VERBOSE,
			"exe=s"					=> \$EXE,
			"bin=s"					=> \$BIN,
			"cleanup"				=> \$CLEANUP,
		  );

$INFILE = shift @ARGV;
$OUTFILE = $INFILE unless ($OUTFILE);
my ($filetag) = $INFILE =~ /(.+)\.gbk/;

# import configuration directly
if ($CONF && -e $CONF) {
	$setup = XMLin($CONF);
}
			
# append underscore to locus tag prefixes ending in a numeric
$LOCUS_TAG_PREFIX = $filetag unless ($LOCUS_TAG_PREFIX);
$LOCUS_TAG_PREFIX = $LOCUS_TAG_PREFIX . '_' if ($LOCUS_TAG_PREFIX =~ /(\d$|X|x)/);
# $LOCUS_TAG_PREFIX =~ s/(\d)$/$1_/; # Does this not work?

unless ($EXE) {
	$EXE = "g3-from-scratch.csh $filetag.fsa $filetag";  #default
#	$EXE = "g3-iterated.csh $filetag.fasta $filetag";  #default
}
else {
	$EXE .= " $filetag.fsa $filetag";
}

if ($BIN) {
	$BIN .= '/' unless ($BIN =~ /$\//);
	$EXE = $BIN . $EXE;
}

unless ($INFILE) {
	die "No file for input specified.\n";
}


# SeqIO input object
my $seqi = Bio::SeqIO->new(
	-file		=> $INFILE,);
my $seq = $seqi->next_seq();

# Use of common temp seq here would be useful
unless (-e "$filetag.fsa") {
	my $tempseq = Bio::SeqIO->new(
		-file		=> ">$filetag" . ".fsa",
		-format		=> "fasta");
	$tempseq->write_seq($seq);
}


my $seqo = Bio::SeqIO->new(
		-format	=> "genbank", 
		-file 	=> ">$OUTFILE");
$seqo->write_seq($seq);

if ($CLEANUP) {
	print "Cleanning up...\n" if ($VERBOSE);
	# Cleanup
	unlink "$INFILE.temp";
	unlink "$filetag.longorfs";
	unlink "$filetag.train";
	unlink "$filetag.icm";
	unlink "$filetag.predict";
	unlink "$filetag.detail";
}



