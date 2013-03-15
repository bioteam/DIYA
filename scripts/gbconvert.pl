#!/usr/bin/env perl

=head1 NAME

gbconvert.pl

=head1 DESCRIPTION

Script to create the tabular and fasta files required by the NCBI
application tbl2asn, then run tbl2asn to create the files 
that NCBI Genome wants.

=head1 USAGE

Example:

/Jake/apps/DIYA/scripts/gbconvert.pl -a WQG -t 47735 -host 'Bos taurus' \
-co 'United States:NE' -cd 10-Nov-2010 -is 'lung' -gc 12 -as 'Celera 7.0' \
-sn 'Isolated from cattle lung, post-mortem' -st '454 Titanium, PacBio RS' \
t/data/2012_03_11_14_27_03-MARC::cmsearch.out.gbk

=cut

use strict;
use diya::MARC::GenbankConvertUtil;
use Getopt::Long;
use Bio::SeqIO;
use FileHandle;
use Bio::Annotation::Collection;
use Bio::Annotation::Comment;

my (
    $debug,            $help,             $project,
    $qual,             $agp,              $taxid,
    $accession_prefix, $host,             $country,
    $collection_date,  $isolation_source, $submission_note,
    $Assembly_Name,    $Sequencing_Technology, $topology
);

# Defaults
my $template        = '/Jake/apps/DIYA/template';
my $executable      = '/Jake/apps/bin/tbl2asn';
my $gcode           = '11';
my $Assembly_Method = '454 Titanium, PacBio RS';
my $Genome_Coverage = '12';

GetOptions(
    "project|p=i" => \$project,
    "help!"       => \$help,
    "qual|q=s"    => \$qual,
    "agp!"        => \$agp,
    "taxid|t=i"   => \$taxid,
    "debug!"      => \$debug,
    "template=s"  => \$template,
    "a=s"         => \$accession_prefix,
    "e=s"         => \$executable,
    "host=s"      => \$host,
    "co=s"        => \$country,         
    "cd=s"        => \$collection_date,
    "is=s"        => \$isolation_source,
    "sn=s"        => \$submission_note, 
    "gc=s"        => \$gcode,
    "an=s"        => \$Assembly_Name,
    "gc=s"        => \$Genome_Coverage, 
    "st=s"        => \$Sequencing_Technology, 
    "as=s"        => \$Assembly_Method,
    "topology=s"  => \$topology
);

usage() if $help;

die "File $template not found" if ( ! -e $template );
die "Executable $executable is not found or not executable" if ( ! -x $executable );
die "Locus tag is required"             if ( ! $accession_prefix );
die "Host is required"                  if ( ! $host );
die "Taxonomy id is required"           if ( ! $taxid );
die "Country is required"               if ( ! $country );
die "Collection date is required"       if ( ! $collection_date );
die "Isolation source is required"      if ( ! $isolation_source );
die "Submission note is required"       if ( ! $submission_note );
die "Sequencing Technology is required" if ( ! $Sequencing_Technology );

my $parser = diya::MARC::GenbankConvertUtil->new(
    -template              => $template,
    -executable            => $executable,
    -accession_prefix      => $accession_prefix,
    -host                  => $host,
    -country               => $country,
    -collection_date       => $collection_date,
    -isolation_source      => $isolation_source,
    -submission_note       => $submission_note,
    -Genome_Coverage       => $Genome_Coverage,
    -Sequencing_Technology => $Sequencing_Technology,
    -Assembly_Method       => $Assembly_Method,
    -Assembly_Name         => $Assembly_Name,
    -gcode                 => $gcode,
    -taxid                 => $taxid,
    -topology              => $topology
);

$parser->debug($debug) if defined $debug;

my $infile = shift @ARGV or usage('Need a Genbank format file');
my $in = Bio::SeqIO->new(-file => $infile,-format => 'genbank');
my $seq = $in->next_seq;

my $id = $seq->id;
$parser->id($id);
my $outdir = $parser->outdir($id);
my $definition = $parser->edit_definition($seq->desc);
$seq->desc($definition);

my $tblfn = "$outdir/$id.tbl";
my $outfeat = FileHandle->new(">$tblfn") or warn ("$tblfn: $!");
my $fsafn = "$outdir/$id.fsa";
my $outfsa = Bio::SeqIO->new(-file => ">$fsafn",
			     -format => 'fasta') or warn ("$fsafn: $!");

$parser->make_namemap($infile);

my @oldFeatures = $seq->remove_SeqFeatures;

# Make fasta and *tbl files
$parser->fixAndPrint(\@oldFeatures,$outfeat,$outfsa,$seq,$definition);

# Calculate overall coverage and add a 'source' feature -
# if the Genbank file is external there's no coverage data
my $comment = $parser->make_top_comment;

if ( $comment ) {
	my $ann = Bio::Annotation::Comment->new;
	$ann->text($comment);
	my $coll = new Bio::Annotation::Collection;
	$coll->add_Annotation('comment',$ann);
	$seq->annotation($coll);
}

# Sort features by location
my @newFeatures = sort { $a->start <=> $b->start } $parser->newFeatures;
$seq->add_SeqFeature(@newFeatures);

$parser->edit_asn_file($seq->desc);

# This first tbl2asn run creates a "discrp" file that has to be parsed
$parser->run_tbl2asn($comment,1);

# Read the discrp file and write a new *.tbl file, fixing discrepancies
$parser->fix_discrp;

# Create a quality file for tbl2asn
$parser->create_qual($qual) if $qual;

$parser->create_agp($infile) if $agp;

$parser->create_cmt();

# Run again
$parser->run_tbl2asn($comment,2);

$parser->cleanup;

sub usage {
	my ($message) = @_;
	$message = "$0 [-h] gbfile" if (! $message );
	print "$message\n";
	exit(-1);
}

__END__

