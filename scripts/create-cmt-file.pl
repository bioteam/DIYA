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

=cut

use strict;
use warnings;
use Getopt::Long;
use LWP::UserAgent;

my ($id,$html);
my $url = 'http://moose/minilims/plugins/MIGS/export_mims.php?instance';

GetOptions( "i=s" => \$id );

die "Need DIYA id" unless ($id);

my $ncbidir = "${id}-gbsubmit.out.d";
my $file    = "${id}.cmt";
my $ua = LWP::UserAgent->new;

my $response = $ua->get("$url=$id");

if ( $response->is_success ) {
    $html = $response->decoded_content;
}
else {
    die $response->status_line;
}

my ($data) = $html =~ /START##\n(.+)\n##MIGS/s;
my (@lines) = $data =~ /(.+)\n/g;
my $text = "StructuredCommentPrefix ##MIGS-Data-START##\n";

for my $line ( @lines ) {
    my ($key,$val) = $line =~ /^(\S+)\s+::\s+(.+)/;
    $val =~ s/_/ /g;
    $text .= "$key\t$val\n";
}

`mkdir $ncbidir` if ( ! -d $ncbidir );
open MYIN,">$ncbidir/$file";
print MYIN $text;
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
