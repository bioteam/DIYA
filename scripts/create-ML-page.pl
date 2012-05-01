#!/usr/bin/perl

=head1 NAME

create-ML-page.pl

=head1 DESCRIPTION

Extracts some information from a DIYA GenBank file and creates a
MiniLIMS page for the DIYA annotation.

=head1 DIYA step and loading script

  <script>
    <name>create-ML-page</name>
    <executable>create-ML-page.pl</executable>
    <command>-i MYSEQID -d OUTPUTDIR -g INPUTFILE.gbk</command>
    <home>/Jake/apps/DIYA/scripts</home>
    <inputformat></inputformat>
    <inputfrom>MARC::cmsearch</inputfrom>
  </script>

Usage  ./DIYAAnnotationResult.php -n <page name> \
                                  -a <assembly page> \
                                  -i <sequence file> \
                                  -l <sequence length> \
                                  -c <number of contigs> \
                                  -o <DIYA results> \
                                  -g <Genbank files>
=cut

use strict;
use warnings;
use Bio::SeqIO;
use Getopt::Long;

my ($id,$outdir,$count,$len, $page, $name, $assembly);
my $loadscript = '/Jake/apps/bin/DIYAAnnotationResult.php';

GetOptions( "i=s" => \$id,
            "d=s" => \$outdir ,
	    "m=s"=>\$page,
	    "s=s"=>\$assembly	    
	    );

die "Need DIYA id, and output directory name" 
  unless ( $outdir && $id );

my $ncbidir = "${id}-gbsubmit.out.d";
 $name="${id}_Annotation";
my $gbk     = "$ncbidir/$id.gbf";

my $in = Bio::SeqIO->new( -file => $gbk,
			  -format => 'genbank' );

while ( my $seq = $in->next_seq ) {
    $count++;
    $len += $seq->length;
}

my $cmd = "$loadscript -g $ncbidir -o $outdir -c $count -l $len -n '$name' -m '$page' -a '$assembly' -i '$outdir/$id.fna'";
print "$0 loading command is $cmd\n";
my $output =`$cmd`;
print $output;

__END__

