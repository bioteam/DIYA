#!/usr/bin/perl 
# $Id: diya-assemble_pseudocontig.pl 348 2009-07-07 15:45:47Z briano $
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

diya-assemble_pseudocontig.pl 

=head1 SYNOPSIS

% diya-assemble_pseudocontig.pl -infile -seqid [-reference] [-spacer] [-species] 
         [-strain] [-gcode] [-newbler] [-readlen] [-trim] [-eliminate] [-clean]

=head1 DESCRIPTION

The basic idea is to take an input file (generally a set of 454 assemblies) and 
then tile against a reference genome and create a pseudocontig in genbank 
format.  It is envisaged as the first step in the Diya pipeline.  This script 
has been designed to work even if a reference genome is not given (there will be 
no tiling, just the same order that the contigs are in).  

The script will run promer and show-tiling (part of the MUMmer package) if the
-reference option is used.

The script creates various output files, but the key output file (with pseudocontigs,
optionally trimmed and cleaned of low-quality contigs) is in Genbank format (*gbk).

There are three ways to input data to the script.

=over 3

=item -infile

The contigs willl be tiled in the same order laid out in the fasta file.

=item  -infile with -reference

promer will be used to order the contigs against a reference sequence.

=item  -infile with -scaff

A *.scaff file will be used to tile the contigs.

=back

If both -reference and -scaff are specified, then -scaff will take precedence.

=head1 PARAMETERS

=over 14

=item --infile

Fasta file, usually from a Newbler assembler run, like 454LargeContigs.fna

=item --reference

The input sequence assemblies will be aligned against this nucleotide
reference genome using promer.

=item --trim

Remove Ns at ends of contigs

=item --eliminate

Remove low quality contigs which have greater than saome percentage of N. Also remove
contigs which are too short. Use this option when using --reference. Supply 2 numbers
separated by a colon, <length>:<percentage cutoff>, like 250:20.

=item --poffset

This is a switch that allows the offset calculated by promer to be 
incorporated into the pseudocontig calculations.  Without this switch, the 
offset will not be recorded.  This option will be meaningless unless the 
-reference option is used.

=item --scaff

A tab-delimited scaffolding file in the format
<contigID>	<direction(EB or BE)> <length of contig> <offset (nt)> <comment>
EB stands for end-beginning and BE is the reverse. If the offset is unknown or unsure 
then a "?" character can be used. Comment is optional.

=item --seqid

Name of the project.

=item --spacer

Sets pseudocontig spacer, default = NNNTTAATTAATTAANNN. If you do not want 
a spacer put a value of "" or 0.  An alternative spacer with starts in six 
frames preceded by S-D sequences: 

NCCATTTTTCCTCCATTTTTCCTCCACTTTTCCTCCNNNTTAATTAATTAANNNGGAGGAAAAATGGAGGAAAAATGGAGGAAAAATGGN

=item --species

The species of the organism

=item --strain

Strain name (not the NS id)

=item --gcode

Genetic code, default = 11

=item --newbler 

Specifies that input file is the output of the newbler assembler 
(e.g. 454LargeContigs.fna) and the coverage can be calculated for each contig 
and added as an annotation.  This option can also be used to 
set exclusion criteria.  If a value in the format <number>:<number> is given 
(e.g. 250:10) all contigs with a length less than or equal will be excluded 
AND all contigs with calculated coverage less than the second number will be 
excluded.

=item --readlen 

The average readlength of the reads in the assembly.  Default value is 105 nt.

=item --clean

Remove all the output fasta files and indices, leave the Genbank file

=back

=head1 TO DO

Talk with Andy on how to best incorporate date, species, gcode, strain fields etc. 
into the Annotation factory to make the pseudocontig. Get the checks for spaces in 
files names sorted. 

Major - implement an XML-based description of the tiling to allow for spliting 
contigs, offsets etc. 

=head1 AUTHOR

Tim Read, timothy.read@med.navy.mil

=cut

use strict;
use Bio::DB::Fasta;
use Getopt::Long;
use Bio::SeqIO;
use Bio::SeqFeature::Generic;

# Variables with default settings
my $SPACER = 'NNNTTAATTAATTAANNN'; # 6-frame translational stop
my $startcoord = 1;
my $READLEN = 105;
my $GCODE = 11;

# Variables that are used to create the Genbank output file 
# Hash of calculated coverage (reads/nucleotide) - the output of the -newbler option block.
my (%copynum,@fastafeats,@cbtfeats);
my $pseudosequence = "";

# Report variables
my ($lengthExcludedcontigs, $total_tiled_contig_length, $total_untiled_contig_length,
	 $No_tiled_ids, $No_untiled, $numExcludedcontigs, $ps_length) = '0';

# Other variables
my (%tiled_against, @tiledids, @remainids, %remain_coverage, %directions, %offset);

# Input variables
my ($INFILE, $OUTDIR, $BIN, $REFERENCE, $POFFSET, $CLEAN,
	 $SEQID, $SPECIES, $STRAIN, $NEWBLER, $SCAFFFILE, $TRIM, $ELIMINATE ) = '';

GetOptions(	"infile=s"				=> \$INFILE,
				"outdir=s"				=> \$OUTDIR,
				"binpath|bin=s"		=> \$BIN,
				"nucmer|reference=s"	=> \$REFERENCE,
				"trim|t"             => \$TRIM,
				"eliminate|e=s"      => \$ELIMINATE,
				"poffset"				=> \$POFFSET,
				"scaff=s"				=> \$SCAFFFILE,
				"seqid|id=s"			=> \$SEQID,
				"spacer=s" 				=> \$SPACER,
				"species=s"				=> \$SPECIES,
				"strain=s"				=> \$STRAIN,
				"gcode=i"				=> \$GCODE,
				"newbler=s" 			=> \$NEWBLER,
				"clean"              => \$CLEAN,
				"readlen=i" 			=> \$READLEN );

# Warnings
check_inputs();

# Fix inputs
$OUTDIR = './' unless $OUTDIR;
$SCAFFFILE = $SEQID . ".scaff" unless $SCAFFFILE;
$BIN .= '/' unless ( $BIN && $BIN =~ /\/$/ );

unless ($SCAFFFILE && -e $SCAFFFILE) {
	print "Scaffold file appears not to exist at $SCAFFFILE. Carrying on without scaffolding\n";
	$SCAFFFILE = '';
}

# The *scaff file supercedes use of the reference, in the event that they are both specified
if ($SCAFFFILE) {
	$REFERENCE = "";
	print "Ignoring -reference option becasue -scaff takes precedence\n"; 
}

# Copy infile to cwd 
system "cp $INFILE $OUTDIR/$SEQID.fna" unless (-e "$OUTDIR/$SEQID.fna");
$INFILE = "$OUTDIR/$SEQID.fna";
chmod 0666, $INFILE;

trim_contigs() if $TRIM;

eliminate_contigs($ELIMINATE) if $ELIMINATE;

my $todays_date = get_date();

my $FASTAFILE = "$OUTDIR/$SEQID.fasta";
my $PSEUDOOUTFILE = "$OUTDIR/$SEQID.gbk";

my $out = Bio::SeqIO->new(-file => ">>$FASTAFILE", -format => 'fasta' );
my $pseudoout = Bio::SeqIO->new(-file => ">>$PSEUDOOUTFILE", -format => 'genbank');

# Create Bio::DB::Fasta database and SeqIO streams, the index file
# is written to the same directory as the fasta file
my $db = Bio::DB::Fasta->new($INFILE);
my @ids = $db->ids;

#-------------------------------------------------------------------------------
# Bio::DB::Fasta does not store the description line so we may need to store it

get_newbler_data() if ( $NEWBLER || $REFERENCE );

#-------------------------------------------------------------------------------
# This section deals with input using either the -reference or -scaff option

if ($REFERENCE) {

	check_for_dna();
	check_for_files();

	my $cmd = $BIN . "promer -p $SEQID $REFERENCE $INFILE"; 
	print "Running $cmd\n";
	system $cmd || die "promer is either not installed or has not been put in the PATH";

	$cmd = $BIN . "show-tiling $SEQID.delta > $OUTDIR/$SEQID.tiling";
	print "Running $cmd\n";
	system $cmd || die 
	  "show-tiling (part of the MUMmer) is either not installed or has not been put in the PATH";

	# Check the byte size of the file to see if its greater than 50 (i.e. just the name).
	# If it fails it's probably becasue the promer output files could not be written.
	if ( (stat("$OUTDIR/$SEQID.tiling"))[7] < 50 ) { 
		print "$OUTDIR/$SEQID.tiling file appears not to have written or is too short to be useful. " .  
		  "Check permissions and the original fasta file. Continuing on as if there is no tiling file.\n";
		$REFERENCE = '';
	}

	else {
		open PROMERTILING, "$OUTDIR/$SEQID.tiling"
		  or die "Cannot open \"$OUTDIR/$SEQID.tiling\" file for some reason"; 
		open PROMERSCAFF, ">$OUTDIR/$SEQID.scaff" 
		  or die "Cannot open \"$OUTDIR/$SEQID.scaff\" file for some reason";

		while (<PROMERTILING>) {
			chomp;
			my @F = split(/\t+/);
			my $referencetext = "";
			my $scaffdir;

			if (/\>/) {  # i.e. the line describes the reference fasta file tiled against
				if ($F[0] =~ /(NC\w+)/) {
					# parse the refseq accession number, if on the description line
					$referencetext = "Promer tiling against accession, $1";
				}  else {
					$referencetext = "Promer tiling"
				}
			}
			else {
				push (@tiledids,$F[7]);
				$directions{$F[7]} = $F[6];
				$tiled_against{$F[7]} = $referencetext;

				if ($POFFSET) {
					$offset{$F[7]} = $F[2];
				}
				else {
					# set the value of the offset to a default from the results of the promer tiling	
					$offset{$F[7]} = "?";
				}

				die "$_\nFormat of all lines of tiling array should be <contig name><tab><+ or ->"
				  unless ($directions{$F[7]}	=~ /\-|\+/);

				if ($F[6] eq "+") {
					$scaffdir = "BE";
				} else {
					$scaffdir = "EB";
				}
				my $contigoff = $offset{$F[7]};
				print PROMERSCAFF "$F[7]\t$scaffdir\t$F[3]\t$offset{$F[7]}\t$referencetext\n";
			}
		}
	}
} # close REFERENCE LOOP

if ( $SCAFFFILE ) {
	open SCAFF, $SCAFFFILE || die "Cannot open $SCAFFFILE file for some reason";
	while (<SCAFF>) {
		my @G = split(/\t+/);
		my $scaffdirA;

		die "Check format of .scaff file: appears to be too few columns" if (@G < 4);

		push (@tiledids,$G[0]);
		if ($G[1] eq "EB") {
			$scaffdirA = "+";
		} elsif ($G[1] eq "BE") {
			$scaffdirA = "-";
		} else { 
			die "The .scaff file does not appear to be in the correct format";
		}
		$directions{$G[0]} = $scaffdirA;	
		$tiled_against{$G[0]} = $G[4];	
		$offset{$G[0]} = $G[3];
	}
}


if ( $REFERENCE || $SCAFFFILE ) {
	print "Creating tiled pseudocontig\n";
	# Check whether all these tiled ids are actually in the fasta file - 
	# if not something is probably wrong with the tiling file
	for my $particle ( @tiledids ) {
		my $flag2 = 0;
		for my $particle2 (@ids) {
			$flag2++ if ($particle eq $particle2);
		}

		if ( $flag2 == 0 ) {
			die "The display_id \"$particle\" from the -tiling file does not match the fasta -infile."
           . "Probably the tiling file is from another project or assembly than the fasta file.";
		}
	}

	# Output to tiled fasta file
	for my $descrip  (@tiledids) {
		my $seq_obj = $db->get_Seq_by_id($descrip);

		if ($directions{$descrip} eq "+") {
			$out->write_seq($seq_obj);
			add_to_pseudo($seq_obj,$descrip,1);
			$total_tiled_contig_length += $seq_obj->length;
		}
		elsif ($directions{$descrip} eq "-") { 
			my	$rev = $seq_obj->revcom;
			add_to_pseudo($rev,$descrip,-1);
			$total_tiled_contig_length += $rev->length;

			my $rev_name  = $rev->display_id . "c";
			$rev->display_id($rev_name); # For some reason I cant assign this id
			$out->write_seq($rev);
		}
	}
}

#-------------------------------------------------------------------------------
# Now deal with the contigs that were not tiled in the block above

for my $element (@ids) {
	my $flag = 0;
	for my $element2 (@tiledids) {
		if ($element eq $element2) {
			$flag++;
		}
	}
	if ($flag == 0) {
		push (@remainids,$element);
		if ($NEWBLER) { # used to sort the remains by coverage
			$remain_coverage{$element} = $copynum{$element};
		}
	}
}

if ($NEWBLER) {
   my ($lencutoff, $covcutoff);

	# Parse out the length and coverage cut-offs
	if ( $NEWBLER =~ /(\d+):(\d+)/ ) {
		($lencutoff,$covcutoff) = ($1,$2);
		print "Newbler contig length cutoff = $lencutoff, coverage cutoff = $covcutoff\n";
	}
	# The idea here is to order the remaining untiled contigs based on coverage, from low to high
	my @sorted_remains = (sort { $remain_coverage{$a} <=> $remain_coverage{$b} } 
						keys %remain_coverage);

	for my $rid (@sorted_remains) {

		my $rseq_obj = $db->get_Seq_by_id($rid);

		if ( $copynum{$rid} < $covcutoff || $rseq_obj->length < $lencutoff ) {

         print "$rid excluded from pseudocontig: Length = " . $rseq_obj->length . 
			  " coverage = ". $copynum{$rid} ."\n";
         $numExcludedcontigs++;
         $lengthExcludedcontigs += $rseq_obj->length;

		} else {

         $out->write_seq($rseq_obj);
         $tiled_against{$rid} = "untiled";
         add_to_pseudo($rseq_obj,$rid,1);
         $total_untiled_contig_length += $rseq_obj->length;
         $No_untiled++;

		}
	}
	# Essentially this block does as above as above but without newbler specified information		
} else {

	for my $rid (@remainids) {
		my $rseq_obj = $db->get_Seq_by_id($rid);
		$out->write_seq($rseq_obj);
		$tiled_against{$rid} = "untiled";
		add_to_pseudo($rseq_obj,$rid,1);
		$total_untiled_contig_length += $rseq_obj->length;
		$No_untiled++;
	}

}

#-------------------------------------------------------------------------------
# Before finishing, chop of the terminal cbt, because we want things to look good

$pseudosequence = substr($pseudosequence,0,-(length($SPACER)));
pop(@cbtfeats);

#-------------------------------------------------------------------------------
# Output pseudocontig - not happy about the fact that annotation features don't work

$ps_length = length($pseudosequence);
my $pseudo_id = "$SEQID";
my $pseudodesc ="[organism=$SPECIES] [strain=$STRAIN] [gcode=$GCODE] [date=$todays_date]";
my $pseudo_obj = Bio::Seq->new(-seq => $pseudosequence,
										 -display_id => $pseudo_id,
										 -desc => $pseudodesc,
										 -accession_number => $SEQID);
$pseudo_obj->add_SeqFeature(@fastafeats);
$pseudo_obj->add_SeqFeature(@cbtfeats);
$pseudoout->write_seq($pseudo_obj);

#-------------------------------------------------------------------------------
# Report 
$No_tiled_ids = @tiledids;
# $No_untiled = @remainids;

print "\nTotal size of pseudocontig = $ps_length\n";
print "Contigs tiled  = $No_tiled_ids\n";
print "Total length of tiled contigs = $total_tiled_contig_length\n";
print "Untiled contigs = $No_untiled\n";
print "Total length of untiled contigs = $total_untiled_contig_length\n";
print "Number of excluded contigs = $numExcludedcontigs\n";
print "Total length of excluded contigs = $lengthExcludedcontigs\n";

#-------------------------------------------------------------------------------
# Clean up

if ( $CLEAN ) {
	unlink "$SEQID.fna.index";
	unlink "$SEQID.fna";
	unlink "$SEQID.fasta";
	unlink "$SEQID.fna.tmp";
}

#-------------------------------------------------------------------------------

sub add_to_pseudo {
	my $sequence = $_[0]->seq();

	if (length($pseudosequence) != 0) {
		$pseudosequence = $pseudosequence.$sequence.$SPACER;	
		calc_fasta_record(length($sequence),length($pseudosequence),$_[0]->display_id,$_[2],
								$tiled_against{$_[1]});
	}

	else {
		$pseudosequence = $sequence.$SPACER;	
		calc_fasta_record(length($sequence),length($pseudosequence),$_[0]->display_id,$_[2],
								$tiled_against{$_[1]});	
	}
}

#-------------------------------------------------------------------------------

sub calc_fasta_record {
	# Stores the information for the fastarecords features

	my $startcoord = 1+$_[1]-$_[0]-(length($SPACER));
	my $endcoord = $_[1] - (length($SPACER));
	if ($SPACER) {
		# create_cbt($startcoord,$startcoord+length($SPACER)); - 
		# this creates the cbt before the added sequence
		create_cbt($_[1]-length($SPACER)+1, $_[1]);	# this creates the cbt after the added sequence
	}
	my $fastarecordfeat = new Bio::SeqFeature::Generic(-start  => $startcoord,
																		-primary => 'fasta_record',
																		-end    => $endcoord,
																		-tag => {name => "$_[2]"},
																		-strand => $_[3]);

	$fastarecordfeat->set_attributes(-tag => {note => "$_[4]"})	if ($_[4]);

	$fastarecordfeat->set_attributes(-tag => {note => 
	  "Calculated sequence coverage = $copynum{$_[2]} reads per base"}) if ($copynum{$_[2]});

	push (@fastafeats, $fastarecordfeat);
}

#-------------------------------------------------------------------------------

sub create_cbt {
	# Creates the contig boundary tag features
	my $cbtrecord = new Bio::SeqFeature::Generic(
    -start  => $_[0],
	 -primary => 'cbt',
    -end    => $_[1],
    -tag => {note => 'contig boundary tag: artificial spacer between sequence contigs',color => 1},
    -strand => 1);

 	push (@cbtfeats, $cbtrecord);
}

#-------------------------------------------------------------------------------

sub get_date {
	my @ds = localtime(time);
	my $date = $ds[4] + 1 . "\-" . $ds[3] . "\-" . ($ds[5] + 1900);
	$date;
}

#-------------------------------------------------------------------------------

sub trim_contigs {
	my $tmpfile = $INFILE . ".tmp";
	`rm $tmpfile` if ( -e $tmpfile );
	`mv $INFILE $tmpfile`;
	my $in = Bio::SeqIO->new(-file => $tmpfile, 
									 -format => 'fasta',
									 -alphabet => 'dna' );
	my $out = Bio::SeqIO->new( -file => ">$INFILE", -format => 'fasta' );

	while ( my $seq = $in->next_seq ) {
		my $str = $seq->seq;
		$str =~ s/^N+//i;
		$str =~ s/N+$//i;
		if ( $str ) {
			$seq->seq($str);
			$out->write_seq($seq);
		}
	}
}

#-------------------------------------------------------------------------------

sub eliminate_contigs {
	# If the percent N of the contig is greater than the 2nd argument or length of 
	# the contig is less than the first argument then eliminate the contig. Example
	# input: 250:10
	my $params = shift;
	my ($len,$cutoff) = $params =~ /(\d+):(\d+)/;
	die "Incorrect argument to --eliminate, length is $len, cutoff is $cutoff" 
	  if ( ! $len || ! $cutoff );

	my $tmpfile = $INFILE . ".tmp";
   `rm $tmpfile` if ( -e $tmpfile );
	`mv $INFILE $tmpfile`;
	my $in = Bio::SeqIO->new( -file => $tmpfile, 
									  -format => 'fasta',
									  -alphabet => 'dna' );
	my $out = Bio::SeqIO->new( -file => ">$INFILE", -format => 'fasta' );

	while ( my $seq = $in->next_seq ) {
		my (@ns) = $seq->seq =~ /(n)/ig;
		if ( ! $seq->length || ( $#ns + 1 )/($seq->length) >= ($cutoff * .01) || $seq->length < $len ) {
			print "Eliminating contig " . $seq->display_id . ", percent N > $cutoff or length < $len\n";
			next;
		}
		$out->write_seq($seq);
	}
}

#-------------------------------------------------------------------------------

sub check_for_dna {
	# Test that reference data is in fasta format and is all DNA
	my $refseqin = Bio::SeqIO->new( -format => 'fasta' , -file => $REFERENCE );

	while( my $refseqobj = $refseqin->next_seq ) { 
		die "The reference file appears to contain sequence that is not DNA" 
		  unless ( $refseqobj->alphabet eq 'dna' );
	}
}

#-------------------------------------------------------------------------------

sub check_for_files {
	# Test for existing output files to avoid accidental overwriting
	die "There is a MUMmer/promer/nucmer out file called $OUTDIR/$SEQID.delta in the directory." 
	  . "Please delete or rename this file before resuming" if (-e "$OUTDIR/$SEQID.delta");

	die "There is a MUMmer/promer/nucmer out file called $OUTDIR/$SEQID.tiling in the directory."
      . "Please delete or rename this file before resuming" if (-e "$OUTDIR/$SEQID.tiling");

	# Test to make sure there are no spaces in the $INFILE or $REFERENCE names,
	# this will trip up promer
	die "The -infile file \"$INFILE\" appears to contain a space character. Please rename the file"
	  . "to contain no spaces and try again" if ($INFILE =~ /\s+/);

	die "The -reference file \"$REFERENCE\" appears to contain a space character.  Please rename"
	  . "the file to contain no spaces and try again"	if ($REFERENCE =~ /\s+/);
}

#-------------------------------------------------------------------------------

sub check_inputs {
	die "Must specify an input file in fasta format using the -infile option at the command line"
	  unless $INFILE;

	die "Must specify a project name using the -seqid option at the command line"
	  unless ($SEQID);

	die "The project name \"$SEQID\" appears to contain a space character. " . 
	  "Please rename the file to contain no spaces and try again" if ($SEQID =~ /\s+/);
}

#-------------------------------------------------------------------------------

sub get_newbler_data {
	open NEWBLERIN,$INFILE;

	while (<NEWBLERIN>) {

		if ( />(\S+).+?length=(\S+)\s+numreads=(\S+)/i ) {
			my ($id,$len,$numreads) = ($1,$2,$3);

			# Skip anything with a length of 0
			next  if ( ! $len );

			my $CN = ($numreads/$len) * $READLEN;
			$copynum{$id} = sprintf("%.2f",$CN);
		}
	}
	die "The fasta file $INFILE is probably not the output from the newbler assembler"
		  . " and the program will not be able to calculate the reads per base."
			 unless ( keys %copynum );
}

#-------------------------------------------------------------------------------

