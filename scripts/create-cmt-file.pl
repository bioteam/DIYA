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
use HTML::Parser;

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




__END__

