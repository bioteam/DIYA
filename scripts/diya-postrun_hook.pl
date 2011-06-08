#!/usr/bin/perl -w
# $Id: diya-postrun_hook.pl 176 2008-08-07 17:34:53Z briano $
#--------------------------------------------------------------------------
# ©Copyright 2008
#
# This file is part of DIYA.
#
# DIYA is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# DIYA is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with the diya package.  If not, see <http://www.gnu.org/licenses/>.
#--------------------------------------------------------------------------

#$ -cwd
#$ -S /usr/bin/perl 
#$ -V 

# This file is intended to serve as the template for a hook to be executed
# upon completion of a diya run.  The project name and path to the output
# file are provided.

use lib '/RemotePerl';
use lib '/RemotePerl/5.8.6/darwin-thread-multi-2level/';

my ($CONF, $SEQID, $FILEPATH, $DIYAPATH) = @ARGV;

my $host = "loki.bdrd";

#-------------------------------------------------------------------------------
#  faa file protein extraction
print "diya-postrun: extracting proteins...\n";
#system $DIYAPATH . "/aux/diya-extract_proteins.pl " . $FILEPATH;

# copy to /shared/data
print "\tcopying to /shared/data...\n";
#`cp $SEQID.faa /shared/data`;
#system "/shared/bin/formatdb -p T -i /shared/data/$SEQID.faa";


#-------------------------------------------------------------------------------
#  xml genome report
#print "diya-postrun: generating genome stat report...\n";
system $DIYAPATH . "/aux/diya-report.pl " . $FILEPATH;

#-------------------------------------------------------------------------------
#  Wikilims-diya.
print "diya-postrun: generating genome wiki page...\n";

use Perlwikipedia;

my $bot_host = $host;
my $bot_user = "WikilimsBot";
my $bot_pass = "WikilimsBot!";

# update genome project page

#	load bot
my $bot = Perlwikipedia->new();
$bot->set_wiki($bot_host,'wiki');

my $loginError = $bot->login($bot_user, $bot_pass);
if ($loginError) {
	error("Login Failure");
	return undef;
}

# deserialize xml stat report
use XML::Simple;
use FileHandle;

my $xmlfn = $FILEPATH;
$xmlfn =~ s/gbk/xml/;


my $xs = XML::Simple->new();
my $xml = $xs->XMLin($xmlfn);

#my $xmlfn = $PROJPATH . '/' . $SEQID . '.xml';
#my $xmlf = FileHandle->new($xmlfn);
#my $xml = XMLin($xmlf);

#my ($fileroot) = $FILEPATH =~ /(\/G_labdata\/454_data\/organized\/annotations\/\d\d\d\d_\d\d\/\w+\/)/;
my $fileroot = `pwd`;
chomp $fileroot;
$fileroot =~ s/$SEQID\.gbk//;
$fileroot =~ s/\/nfs\/Asgard//;
my $urlroot = "http://" . $host . $fileroot;

# create genome page
my $genome_text = "{{Genome\n";	
$genome_text .= "|id="	.						$xml->{genome}->{seqid} . "\n";
$genome_text .= "|species=" . 					$xml->{genome}->{species} . "\n";
$genome_text .= "|strain=" . 					$xml->{genome}->{strain} . "\n";
$genome_text .= "|date=" . 						$xml->{genome}->{date} . "\n";
$genome_text .= "|ns=" . 						$xml->{genome}->{ns} . "\n";
$genome_text .= "|topology=" .					$xml->{genome}->{sequence}->{topology} . "\n";
$genome_text .= "|molecule=" .					$xml->{genome}->{sequence}->{molecule} . "\n";
$genome_text .= "|length=" .					$xml->{genome}->{sequence}->{length} . "\n";
$genome_text .= "|pct_gc=" . 					$xml->{genome}->{sequence}->{gc} . "\n";
$genome_text .= "|pct_coding=" .				$xml->{genome}->{sequence}->{coding} . "\n";
$genome_text .= "|contigs=" . 					$xml->{genome}->{assembly}->{contigs} . "\n";
$genome_text .= "|genes=" . 					$xml->{genome}->{features}->{genes} . "\n";
$genome_text .= "|ncrnas=" . 					$xml->{genome}->{features}->{ncRNA} . "\n";
$genome_text .= "|avg_gene_length=" .			$xml->{genome}->{features}->{avg_gene_length} . "\n";
$genome_text .= "|avg_intergenic_length=" .		$xml->{genome}->{features}->{avg_intergenic_length} . "\n";
$genome_text .= "|overlaps=" . 					$xml->{genome}->{features}->{overlaps} . "\n";

$genome_text .= "|filepath=" .					$FILEPATH . "\n";
$genome_text .= "|fileroot=" .					$fileroot . "\n";
$genome_text .= "|urlroot=" . 					$urlroot . "\n";

$genome_text .= "}}\n";

my $comment = "Reporting genome results from DIYA run";
	
eval {
	$bot->edit($SEQID, $genome_text, $comment) or print "error: $bot->{errstr}\n";
};
if ($@) {
	print "$@" . "\t $bot->{errstr}\n";
}

#		charts?


#-------------------------------------------------------------------------------
# gbrowse
my $gbrowse_host = $host;
my $gbrowse_user = "www";
my $gbrowse_pass = "www";
my $gbrowse_dsn = "\"dbi:mysql:TEST;host=$gbrowse_host\"";

use Bio::DB::GFF;
$db = Bio::DB::GFF->new(
	-adaptor=>'dbi::mysqlopt', 
	-dsn=>'dbi:mysql:TEST;host=loki.bdrd', 
	-user=>'www',
	-pass=>'www',
);

#delete old
eval {
    $db->delete($SEQID);
};
if ($@) {
    use Carp;
    use Data::Dumper;
    Carp::cluck($@);
    print Dumper $db, $SEQID;
}

my $make_gff_cmd = "$DIYAPATH" . "/aux/diya-genbank2gff3.pl";
$make_gff_cmd .= " --seqid $SEQID";
$make_gff_cmd .= " $FILEPATH";
print "diya-postrun: generating gff3 file...\n";
system $make_gff_cmd;

#add new
my $load_gff_cmd = "bp_load_gff.pl";
$load_gff_cmd .= " --user $gbrowse_user";
$load_gff_cmd .= " --pass $gbrowse_pass";
$load_gff_cmd .= " --dsn $gbrowse_dsn";
$load_gff_cmd .= " $FILEPATH.gff";
#
print "diya-postrun: loading Bio::DB::GFF database...\n";
system $load_gff_cmd;
