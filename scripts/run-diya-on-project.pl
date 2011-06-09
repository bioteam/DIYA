#!/usr/bin/perl
# $Id: run-diya-on-project.pl 272 2010-01-12 22:22:28Z briano $

=head1 NAME

run-diya-on-project

=head1 DESCRIPTION

Script to run the DIYA pipeline on BDRD genomes. Example:

 ./run-diya-on-project.pl -p PMP4 -c genome-w-ref-wo-qual.conf

Or, to do just a subset of strains:

 ./run-diya-on-project.pl -p PMP4 -s NS4522 NS2311 NS3600

=head1 REQUIREMENTS

=over 3

=item A *conf file

=item A genome reference file (fasta format)

=item An ASN "template" file

=back

=cut

use strict;
use MediaWiki::Bot;
use Getopt::Long;
use HTML::Strip;

# Defaults
my $conf = 'genome-w-ref-wo-qual.conf';
# my $reference = '/site/data/reference/NC_005945.fa';
my $reference = 'AE017334.fa';
my $filedir = '/massive/data/genome/454_data/organized/assemblies/';

# The 'wikilims' directory below is the SVN wikilims repository -
# You can check it out like this: svn co file:///massive/site/svn/wikilims
my $nsrowscript = '/site/home/briano/wikilims/trunk/extensions/wikilims/nsrow.pl';

my $project;
my @nsnames;

GetOptions( 'p|project=s'   => \$project,
            'c|conf=s'      => \$conf,
			   'r|reference=s' => \$reference,
			   's=s{1,}'       => \@nsnames );

die "Must supply project name with \'-p\'" if ( ! $project );
die "Reference genome $reference not found" if ( ! -e $reference );
die "$nsrowscript script not found" if ( ! -e $nsrowscript );
die "Directory $filedir not found" if ( ! -d $filedir );
die "Conf file $conf not found" if ( ! -e $conf );

my $wiki = { 'host' => 'loki.bdrd',
             'dir'  => 'wiki',
             'user' => 'WikilimsBot',
             'pass' => 'WikilimsBot!' };

my $bot = MediaWiki::Bot->new;

$bot->set_wiki($wiki->{host}, $wiki->{dir});
$bot->login($wiki->{user}, $wiki->{pass});

my $projectData = get_project_data($project);

for my $genome ( sort keys %{$projectData} ) {

	# skip genomes not sequenced here
	next if $genome !~ /^[a-zA-Z]{5}\d{4}$/;

   # 'bthur0014' => '[[NS5497]] ibl 422 [[Special:Whatlinkshere/NS5497|rev]] [[bthur0014]] [http://www.ncbi.nlm.nih.gov/sites/entrez?db=genomeprj&cmd=Retrieve&dopt=Overview&list_uids=29735  29735] [http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=527020  527020] [[BTH]] Bacillus thuringiensis 1 rdirs 1 pdirs [[P_2008_06_15_14_10_25_loki]] 131911276 bases 523730 reads Bacillus thuringiensis IBL4222 2008_06 '

   my ($ns) = $projectData->{$genome} =~ /\[\[(NS\d+)\]\]\s+/
     or die "No NS name for $genome";

   next if ( @nsnames && ! grep /^$ns$/,@nsnames );

   my ($strain) = $projectData->{$genome} =~ /\[\[NS\d+\]\]\s+(.+?)\s+\[\[/
     or die "No strain name for $genome";
   my ($pdir) = $projectData->{$genome} =~ /pdirs\s+\[\[([^]]+)/
     or die "No P dir for $genome";
   my ($tdir) = $projectData->{$genome} =~ /(\d{4}_\d{2})\s*$/
     or die "No T dir for $genome";
   my ($tid) = $projectData->{$genome} =~ /wwwtax.cgi\?id=(\d+)/
     or die "No NCBI taxon id for $genome";
   my ($gspecies) = $projectData->{$genome} =~ /\[\[[A-Z]{3}\]\]\s+(.+?)\s+\d+\s+rdirs/
     or die "No species for $genome";

	my $cmd =
"nohup /site/perl/diya.pl -conf $conf -set MYREF=$reference -set MYUNIREFD=/site/data/uniref50.fasta -set MYTAXID=$tid -set MYSTRAIN=\'$strain\' -set MYSPECIES=\'$gspecies\' -set MYRPSD=/site/data/Cdd -set MYSEQID=$genome -set MYCLUSTERS=/site/data/Clusters.bcp -set MYCDD=/site/data/cddid_all.tbl -set MYPROJECT=$project $filedir/$tdir/$pdir/454AllContigs.fna > $genome.out &";

	`$cmd`;

}


sub get_project_data {
   my $project = shift;
   my $projectData;

   my $text = $bot->get_text($project);

   my ($nstring) = $text =~ /\|\s*NSIDs\s*=\s*([^}|]+)/s;
   my (@nsids) = $nstring =~ /(NS\d+)/g;
   my $nslist = join ',',@nsids;

   my $projectHTML = `$nsrowscript $nslist`;

   my $hs = HTML::Strip->new();
   my $cleanText = $hs->parse($projectHTML);
   $hs->eof;

   # positive lookahead...
   my @lines = split /(?=\[\[NS\d+\]\].+?)/, $cleanText;

   for ( @lines ) {
      /\[\[([a-z]{5}\d{4})\]\]/i;
      $projectData->{$1} = $_;
   }
   $projectData;
}

__END__

nohup /site/perl/diya.pl -conf genome-w-ref-wo-qual.conf -set MYREF=AE017334.fa \
-set MYSTRAIN='m1293' -set MYTAXID=526973 -set MYSPECIES='Bacillus cereus' -set MYSEQID=bcere0001 \
-set MYCLUSTERS=/site/data/Clusters.bcp -set MYCDD=/site/data/cddid_all.tbl -set MYPROJECT=PMP4 \
/Asgard/G_labdata/454_data/organized/assemblies/2008_03/P_2008_03_09_15_37_45_loki/454AllContigs.fna \
> bcere0001.out &

Example P dir:

/massive/data/genome/454_data/organized/assemblies/2007_10/P_2007_10_20_04_31_54_loki/


