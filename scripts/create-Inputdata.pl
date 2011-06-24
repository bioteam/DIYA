#!/usr/bin/perl

=head1 NAME

create-Inputdata.pl

=head1 DESCRIPTION

Script to create the Input.data file for sRNAscanner. Also copies the input sequence
file to the sRNAscanner dir and creates the *ptt file for sRNAscanner.

=cut

use strict;
use Bio::SeqIO;
use Bio::FeatureIO;
use File::Spec;

my $sRNAscannerdir = '/usr/local/share/apps/sRNAscanner';
my $inputdatafile  = 'Input.data';

my $id        = shift or die "Must provide a sequence id";
my $gbkinfile = shift or die "Must provide a Genbank file";
my $outputdir = shift or die "Must provide current working directory";

$outputdir = File::Spec->rel2abs($outputdir);
print "Absolute path of OUTPUTDIR is $outputdir\n";

# clean out the working directories
`rm -fr $sRNAscannerdir/output_files/*`;
`rm -fr $sRNAscannerdir/temp/*`;

my $pttfile = makePtt($id,$gbkinfile);

makeInputdata($id,$pttfile);



sub makePtt {
	my ($id,$gbkinfile) = @_;
	my @fs;

	my $in = Bio::SeqIO->new(-file => $gbkinfile, -format => 'genbank');
	# Get features from the preceding step
	while (my $s = $in->next_seq) {
	 	for my $f ($s->get_SeqFeatures) {
			push @fs, $f;
		}
	}

	# Write features to a *ptt file in the OUTPUTDIR
	my $file = "$outputdir/$id.ptt";
	my $fout = Bio::FeatureIO->new(-file => ">$file", -format => 'ptt');
		for my $f (@fs) {
   		$fout->write_feature($f);
 	}
 	print "*ptt file is $file\n";
 	$file;
}

sub makeInputdata {
	my ($id,$pttfile) = @_;
	my $subfasta = '<FASTA_FILE>';
	my $subptt   = '<PTT_FILE>';
	my $fastain  = "$id.fa";

	$fastain = File::Spec->rel2abs($fastain);
	print "Absolute path of the input fasta is $fastain\n";

	$fastain =~ s|/|\\/|g;
	$pttfile =~ s|/|\\/|g;
	print "Substituted path is $fastain\n";
	print "Substituted path is $pttfile\n";

	`rm $sRNAscannerdir/$inputdatafile` if ( -e "$sRNAscannerdir/$inputdatafile");
	`cp $sRNAscannerdir/$inputdatafile.template $sRNAscannerdir/$inputdatafile`;

	# Substitute strings in Input.data
	`perl -pi.bak -e 's/$subfasta/$fastain/g' $sRNAscannerdir/$inputdatafile`;
	`perl -pi.bak -e 's/$subptt/$pttfile/g' $sRNAscannerdir/$inputdatafile`;
	`rm $sRNAscannerdir/*bak`;
	print "Created $sRNAscannerdir/$inputdatafile\n";
}



__END__

*ptt example:

Salmonella typhimurium LT2, complete genome - 1..4857432
4423 proteins
Location	Strand	Length	PID	Gene	Synonym	Code	COG	Product
190..255	+	21	16763391	thrL	STM0001	-	-	thr operon leader peptide
337..2799	+	820	16763392	thrA	STM0002	-	COG0527E,COG0460E	bifunctional 
aspartokinase I/homeserine dehydrogenase I
2801..3730	+	309	16763393	thrB	STM0003	-	COG0083E	homoserine kinase
