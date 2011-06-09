#!/usr/bin/perl
# $Id: run-diya-on-wiki-genome.pl 192 2009-01-16 14:37:11Z briano $

=head1 NAME

run-diya-on-genome - NAME

gbsubmit.pl

=head1 DESCRIPTION

Script to run the DIYA pipeline on BDRD genomes.

=cut

use strict;
use Perlwikipedia;

my $idfile = shift or die "Must provide a file with genome ids, one per line";

my $conf = 'genome-annotator.conf';

my $fileserver = 'G_labdata/454_data/organized/assemblies/';

my $wiki = { 'host' => 'loki',
				 'dir'  => 'wiki',
				 'user' => 'Bioteam',
				 'pass' => 'bioteampw' };

my $bot = Perlwikipedia->new;
$bot->set_wiki($wiki->{host}, $wiki->{dir});
$bot->login($wiki->{user}, $wiki->{pass});

open MYIN,$idfile or die "Cannot open file $idfile";

while ( my $id = <MYIN>) {
	chomp $id;
	$id = ucfirst($id);

	my $text = $bot->get_text($id);
	my ($pdir) = $text =~ /\|\s*P_dir\s*=\s*(\S+)/;
   my ($strain) = $text =~ /\|\s*strain\s*=\s*(\S+)/;
   my ($species) = $text =~ /\|\s*species\s*=\s*([A-Za-z\s]+)/;
	$strain =~ s/_/ /g;

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
