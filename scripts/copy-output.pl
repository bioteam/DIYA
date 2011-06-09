#!/usr/bin/perl -w
# $Id: copy-output.pl 257 2009-05-27 18:11:12Z briano $
# Script to copy DIYA output to a repository

use strict;
use File::Find;
use File::Basename;

my $dest = '/massive/data/genome/submissions';

my $dir = shift or die "No directory given";

my @predict = <$dir/*predict>;

my ($id) = $predict[0] =~ /([\d\w]+)\.predict$/i or die "No id found";

#my $newid = $id . "0001" if ( $id !~ /\d+$/ );

die "No destination directory" unless ( -d "$dest/$id" );

my $cmd = "cp -f $dir/*gbk $dir/*out $dir/$id* $dir/454AllContigs.fna $dest/$id";

system $cmd;

my @quals;

find( sub { push @quals, $File::Find::name if -f && /$id.qual/ }, "." );

# find( sub { push @quals, $File::Find::name if -f && /$newid.qual/ }, "." );

if ( @quals ){
	$cmd = "cp " . $quals[0] . " $dest/$id";
	system $cmd;
}

if ( -d "$id-gbsubmit.out.d" ) {
	$cmd = "cp $id-gbsubmit.out.d/$id.* $dest/$id";
	system $cmd;
	$cmd = "cp $id-gbsubmit.out.d/discrp $dest/$id";
	system $cmd;
}



__END__

loki 109 ~/genome-annotation/done>ls yaldo-diya/
  2008_11_21_15_15_21-glimmer3.out.gbk 454AllContigs.fna
  2008_11_21_15_16_31-extractCDS.out yaldo.detail
  2008_11_21_15_16_56-blastpCDS.out yaldo.fa
  2008_11_21_15_16_56-blastpCDS.out.gbk yaldo.fasta
  2008_11_21_15_16_56-blastpCDS.out.index yaldo.fna
  2008_11_21_15_47_20-rpsblastCDS.out yaldo.fna.index
  2008_11_21_15_47_20-rpsblastCDS.out.gbk yaldo.gbk
  2008_11_21_15_47_20-rpsblastCDS.out.index yaldo.icm
  2008_11_21_15_53_20-tRNAscanSE.out yaldo.longorfs
  2008_11_21_15_53_20-tRNAscanSE.out.gbk yaldo.out
  2008_11_21_15_54_17-rnammer.out yaldo.predict
  2008_11_21_15_54_17-rnammer.out.gbk yaldo.train

