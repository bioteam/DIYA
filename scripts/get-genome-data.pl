#!/usr/bin/perl
# $Id: get-genome-data.pl 192 2009-01-16 14:37:11Z briano $

=head1 NAME

get-genome-data.pl

=head1 DESCRIPTION

Script to collect genome data.

=cut

use strict;
use Perlwikipedia;

my $idfile = shift or die "Must provide a file with genome data, one line per genome";

my $wiki = { 'host' => 'loki',
				 'dir'  => 'wiki',
				 'user' => 'Bioteam',
				 'pass' => 'bioteampw' };

my $bot = Perlwikipedia->new;
$bot->set_wiki($wiki->{host}, $wiki->{dir});
$bot->login($wiki->{user}, $wiki->{pass});

open MYIN,$idfile or die "Cannot open file $idfile";

while ( <MYIN>) {
	my @data = split "\t";
	# $id = ucfirst($id);

	my ($ns,$strain,$id,$pdir) = ($data[0],$data[1],$data[3],$data[8]);
	$strain =~ s/_/ /g;

	next if ( $pdir !~ /^P/ );

	# need "taxid","genusspecies","ncbiprojectid","species" from NS page
	my $txt = $bot->get_text($ns);
	my ($taxid) = $txt =~ /\|\s*taxid\s*=\s*(\S+)/;
   my ($species) = $txt =~ /\|\s*genusspecies\s*=\s*([A-Za-z]+\s+[A-Za-z]+)/;
   my ($project) = $txt =~ /\|\s*ncbiprojectid\s*=\s*(\d+)/;
   my ($speciesid) = $txt =~ /\|\s*species\s*=\s*([A-Z]+)/;

	$txt = $bot->get_text($pdir);
	my ($year) = $txt =~ /\|\s*year\s*=\s*(\S+)/;
   my ($month) = $txt =~ /\|\s*month\s*=\s*(\S+)/;
	my $tdir = $year . "_" . $month;
	my ($readlen) = $txt =~ /\|\s*readlen_value\s*=\s*(\d+)/;

	print (join ",", ($id,$ns,$strain,$tdir,$pdir,$taxid,$species,$project,$speciesid,$readlen) );
	print "\n";
}


__END__

Example P dir:
http://localhost:16080/G_labdata/454_data/organized/assemblies/2006_10/P_2006_10_06_10_01_01_runAssembly

{{Genome
|taxid=634
|P_dir=2006_10/P_2006_10_06_10_01_01_runAssembly
|ncbiprojectid=16104
|NSnumber=NS2459
|id=yberc0001
|species=Yersinia bercovieri
|strain=ATCC_43970
|date=7-3-2008
|ns=yberc0001
|topology=linear
|molecule=dna
|length=4328018
|pct_gc=48
|pct_coding=85
|contigs=177
|genes=4092
|ncrnas=65
|avg_gene_length=904
|avg_intergenic_length=184.21
|overlaps=643
|filepath=yberc0001.gbk
|fileroot=/G_labdata/454_data/organized/annotations/2008_07/A_PGL1_2008_07_04_14_36_37
|urlroot=http://loki.bdrd/G_labdata/454_data/organized/annotations/2008_07/A_PGL1_2008_07_04_14_36_37
}}
