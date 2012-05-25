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

my ($data) = $html =~ /(START##.+##MIGS)/;

print "Done\n";


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

