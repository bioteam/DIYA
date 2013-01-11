#!/usr/bin/perl

=head1 NAME

create-cmt-file.pl

=head1 DESCRIPTION

Extracts some MIGS data from a MiniLIMS page and creates a *cmt file for tbl2asn
using that information.

=head1 DIYA step

  <script>
    <name>create-cmt-file</name>
    <executable>create-cmt-file.pl</executable>
    <command>-i MYSEQID</command>
    <home>/Jake/apps/DIYA/scripts</home>
    <inputformat></inputformat>
    <inputfrom></inputfrom>
  </script>

=head2 Example URL

http://moose/minilims/plugins/MIGS/export_mims.php?instance=Bibersteinia_trehalsoi_192_lung_swab

=head2 Example command

./create-cmt-file.pl -i test -s Bibersteinia_trehalsoi_192_lung_swab

=cut

use strict;
use warnings;
use Getopt::Long;
use LWP::UserAgent;
use Cwd 'abs_path';

my ( $id, $html, $strain );
my $url = 'http://moose/minilims/plugins/MIGS/export_mims.php?instance';
# Allowed values in MIGS 3.0
my @investigationTypes =
  qw(eukaryote bacteria_archaea plasmid virus organelle metagenome miens-survey miens-culture);

GetOptions( "i=s" => \$id, "s=s" => \$strain );

die "Need DIYA id and strain" if ( ! $id || ! $strain );

$strain     =~ s/\s+/_/g;
my $ncbidir = "${id}-gbsubmit.out.d";
$ncbidir    = abs_path($ncbidir);
my $file    = "${id}.cmt";
my $ua = LWP::UserAgent->new;

my $response = $ua->get("$url=$strain");

if ( $response->is_success ) {
    $html = $response->decoded_content;
} else {
    die $response->status_line;
}

my ($data) = $html =~ /START##\n(.+)\n##MIGS/s;
my (@lines) = $data =~ /(.+)\n/g;
my $text = "StructuredCommentPrefix\t" . '##MIGS:3.0-Data-START##' . "\n";

for my $line ( @lines ) {
    my ($key,$val) = $line =~ /^(\S+)\s+::\s+(.*)/;
    $val =~ s/_/ /g;
    $val =~ s/\s+/_/ if ( $key eq 'investigation_type' );
    $text .= "$key\t$val\n";
}

$text .= "StructuredCommentSuffix\t" . '##MIGS:3.0-Data-END##' . "\n";
print "Text is:\n$text\n";

`mkdir $ncbidir` if ( ! -d $ncbidir );
open MYIN,">$ncbidir/$file";
print MYIN $text;
print "*cmt file is $ncbidir/$file\n";
print "Problem writing to $ncbidir/$file\n" if ( ! -e "$ncbidir/$file" );
print "Done with $0\n";

__END__

<pre>
##MIGS-Data-START##
biome               ::  organ
biotic_relationship ::  commensal
collection_date     ::  2010-10-10
estimated_size      ::  2500000
geo_loc_name        ::  United_States:Nebraska
health_disease_stat ::  dead
investigation_type  ::  bacteria_archaea
material            ::  mucus
pathogenicity       ::  animal
ploidy              ::  haploid
samp_collect_device ::  swab
specific_host       ::  9913
submitted_to_insdc  ::  N
##MIGS-Data-END##
</pre>

Example *cmt file:

StructuredCommentPrefix ##MIGS-Data-START##
investigation_type bacteria_archaea
project_name Bibersteinia trehalosi 192
collection_date 2008-08-09
lat_lon 35.64N 56E
depth 10cm
alt_elev 100m
country     Canada
environment    Soil
num_replicons  4
ref_biomaterial Not available
biotic_relationship     Free living
trophic_level  autotroph
rel_to_oxygen  Facultative
isol_growth_condt      Not available
sequencing_meth pyrosequencing
assembly Newbler v. 2.3 (pre-release)
finishing strategy draft, 15x, 26 contigs
StructuredCommentPrefix ##Genome-Assembly-Data-START##
Assembly Method Celera assembler v. 6.1
Assembly Name CA_gls454027
Genome Coverage 16.3x
Sequencing Technology 454 Paired-End Titanium
StructuredCommentSuffix ##Genome-Assembly-Data-END##
