#!/usr/bin/perl
# $Id: run-diya-on-genome.pl 224 2009-03-11 13:46:44Z briano $

=head1 NAME

run-diya-on-genome

=head1 DESCRIPTION

Script to run the DIYA pipeline on BDRD genomes. Example:

 ./run-diya-on-genome.pl PMP4-noX-data.txt

=cut

use strict;

# example file: PMP4-noX.csv
my $idfile = shift or die "Must provide a file with genome data, one per line";

my $conf = 'genome-w-ref-wo-qual.conf';

my $server = '/Asgard/G_labdata/454_data/organized/assemblies';

open MYIN,$idfile or die "Cannot open file $idfile";

while ( <MYIN> ) {
	next if /^#/;

	my ($ns,$strain,$rev,$id,$pid,$nid,$shortSpecies,$gspecies,$numRdirs,$numPdirs,$pdir,
		             $numBases,$numReads,$gbtaxname,$tdir) = split ",";
	chomp $tdir;

	next unless ( $id && $strain && $gspecies );

	my $cmd = 
"/site/perl/diya.pl -conf $conf -set MYREF=AE017334.fa -set MYUNIREFD=/site/data/uniref50.fasta -set MYTAXID=$nid -set MYSTRAIN=\'$strain\' -set MYSPECIES=\'$gbtaxname\' -set MYRPSD=/site/data/Cdd -set MYSEQID=$id -set MYCLUSTERS=/site/data/Clusters.bcp -set MYCDD=/site/data/cddid_all.tbl -set MYPROJECT=PMP4 $server/$tdir/$pdir/454AllContigs.fna > $id.out";

	system($cmd);
}


__END__
NS5536,BDRD-ST26,rev,bcere0013,29661,526975,BCE,Bacillus cereus,3 rdirs,2 pdirs,P_2008_02_14_13_25_08_loki,132343857 bases,550158 reads,Bacillus cereus BDRD-ST26,2008_02

nohup /site/perl/diya.pl -conf genome-w-ref-wo-qual.conf -set MYREF=AE017334.fa \
-set MYSTRAIN='m1293' -set MYTAXID=526973 -set MYSPECIES='Bacillus cereus' -set MYSEQID=bcere0001 \
-set MYCLUSTERS=/site/data/Clusters.bcp -set MYCDD=/site/data/cddid_all.tbl -set MYPROJECT=PMP4 \
/Asgard/G_labdata/454_data/organized/assemblies/2008_03/P_2008_03_09_15_37_45_loki/454AllContigs.fna \
> bcere0001.out &

Example P dir:

/Asgard/G_labdata/454_data/organized/assemblies/2006_10/P_2006_10_06_10_01_01_runAssembly
