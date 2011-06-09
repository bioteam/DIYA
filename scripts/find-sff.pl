#!/usr/bin/perl
# $Id: create-genome-pages.pl 273 2010-01-27 21:27:58Z briano $

use strict;
use MediaWiki::Bot;
use File::Find;

my $filedir = '/massive/data/genome/454_data/organized/assemblies';

my $wiki = { 'host' => 'loki.bdrd',
				 'dir'  => 'wiki',
				 'user' => 'WikilimsBot',
				 'pass' => 'WikilimsBot!' };

my $bot = MediaWiki::Bot->new;
$bot->set_wiki($wiki->{host}, $wiki->{dir});
$bot->login($wiki->{user}, $wiki->{pass});

my @pages = $bot->get_pages_in_category('Category:Is_a_genome');

for my $page ( @pages ) {
	my @files;
	my $text = $bot->get_text($page);
	my ($dir) = $text =~ /P_dir\s*=\s*(P.+)/;
	my ($y,$m) = $dir =~ /P_(\d+)_(\d+)/;
	# find( sub { push @files, $File::Find::name if -f && /\.sff$/ }, "$filedir/${y}_${m}/$dir" );
	find( sub { push @files, $File::Find::name if -f && /\.ace$/ }, "$filedir/${y}_${m}/$dir" );
	print "$page - @files\n" if @files;
}


__END__

Example P dir:
http://localhost:16080/G_labdata/454_data/organized/assemblies/2006_10/P_2006_10_06_10_01_01_runAssembly

{{Genome
|SubmittedToNCBI=yes
|FullyAnnotated=yes
|Project_Name=PGL1
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


