#!/sw/bin/perl
# $Id: create-genome-pages.pl 273 2010-01-27 21:27:58Z briano $

=head1 NAME

create-genome-pages.pl

=head1 DESCRIPTION

Script to create Category:Genome pages. Requires a project name like
"PMP4" and a genome name.

This script executes the script nsrow.pl to get Project information.

=head1 USAGE

  ./create-genome-pages.pl -p <project name> -g <genome>

Example:

  ./create-genome.pl -p PMP4 -g bcere0016 -s y -a y

=cut

use strict;
use MediaWiki::Bot;
use Date::Format;
use Bio::SeqIO;
use Getopt::Long;
use HTML::Strip;

my $projectName;
my @genomes;
# defaults
my $submittedToNCBI = 'yes';
my $fullyAnnotated = 'yes';

GetOptions( 'p|project:s'    => \$projectName,
				'genome=s{,}' => \@genomes,
			   's|submitted:s' => sub {
					my ($name,$val) = $@;
					$submittedToNCBI = 'yes' if ( $val =~ /^y/i );
               $submittedToNCBI = 'no' if ( $val =~ /^n/i );
				},
			   'a|annotated:s' => sub {
               my ($name,$val) = $@;
					$fullyAnnotated = 'yes' if ( $val =~ /^y/i );
               $fullyAnnotated = 'no' if ( $val =~ /^n/i );
				} );

my $date = time2str("%e-%L-%Y",time);

my $nsrowscript = '/site/home/briano/wikilims/trunk/extensions/wikilims/nsrow.pl';
die "$nsrowscript script not found" if ( ! -e $nsrowscript );

my $filedir = '/massive/data/genome/submissions';
die "Directory $filedir not found" if ( ! -e $filedir );

die "Must supply project name with \'-p\' and genome with \'-g\'"
  if ( ! $projectName || ! $genomes[0] );

my $wiki = { 'host' => 'loki.bdrd',
				 'dir'  => 'wiki',
				 'user' => 'WikilimsBot',
				 'pass' => 'WikilimsBot!' };

my $bot = MediaWiki::Bot->new;
$bot->set_wiki($wiki->{host}, $wiki->{dir});
$bot->login($wiki->{user}, $wiki->{pass});

my $projectData = get_project_data($projectName);

for my $genome (@genomes) {

	die "Could not find genome $genome in project $projectName" 
	  unless ( defined $projectData->{$genome} );

	# 'bcere0016' => '[[NS2805]] 95/8201 [[Special:Whatlinkshere/NS2805|rev]] [[bcere0016]]
	# [http://www.ncbi.nlm.nih.gov/sites/entrez?db=genomeprj&cmd=Retrieve&dopt=Overview&list_uids=29669  29669]
	# [http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=526979  526979] [[BCE]] Bacillus cereus 4 rdirs
	# 5 pdirs [[P_2008_02_13_17_56_34_loki]] 121769577 bases 487766 reads '

	my ($ns) = $projectData->{$genome} =~ /\[\[(NS\d+)\]\]/
	  or print "No NS strain for $genome\n";
	my ($strain) = $projectData->{$genome} =~ /\[\[NS\d+\]\]\s+(.+?)\s+\[\[/
	  or print "No strain name for $genome\n";
	my ($pdir) = $projectData->{$genome} =~ /pdirs\s+\[\[([^]]+)/
	  or print "No P dir for $genome\n";
	my ($tid) = $projectData->{$genome} =~ /wwwtax.cgi\?id=(\d+)/
	  or print "No NCBI taxon id for $genome\n";
	my ($gspecies) = $projectData->{$genome} =~ /\[\[[A-Z]{3}\]\]\s+(.+?)\s+\d+\s+rdirs/
	  or print "No species for $genome\n";
	my ($pid) = $projectData->{$genome} =~ /list_uids=(\d+)/
	  or print "No NCBI project for $genome\n";

	#my $gbk = "$filedir/$genome/$genome.gbf" or die "No Genbank file found";
	my $gbk = "${genome}-gbsubmit.out.d/$genome.gbf" or die "No Genbank file found";

	my ($pct_gc,$pct_coding,$contignum,$genes,$ncrnas,$avg_gene_len,
		 $avg_intergenic_len,$overlaps,$length,$top,$alphabet) = get_seq_data($gbk);

	my $text ="
{{Genome
|Submitted to SRA=
|SubmittedToNCBI=$submittedToNCBI
|FullyAnnotated=$fullyAnnotated
|Project_Name=$projectName
|taxid=$tid
|P_dir=$pdir
|ncbiprojectid=$pid
|NSnumber=$ns
|id=$genome
|species=$gspecies
|strain=$strain
|date=$date
|ns=$genome
|topology=$top
|molecule=dna
|length=$length
|pct_gc=$pct_gc
|pct_coding=$pct_coding
|contigs=$contignum
|genes=$genes
|ncrnas=$ncrnas
|avg_gene_length=$avg_gene_len
|avg_intergenic_length=$avg_intergenic_len
|overlaps=$overlaps
|filepath=
|fileroot=
|urlroot=
}}";

   $bot->edit($genome,$text,"Create Genome page for $genome");
   print "Created Genome page for $genome\n";
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
	my $cleanText = $hs->parse( $projectHTML );
	$hs->eof;

	# positive lookahead...
	my @lines = split /(?=\[\[NS\d+\]\].+?)/, $cleanText;

	for ( @lines ) {
		/\[\[([a-z]{5}\d{4})\]\]/i;
		$projectData->{$1} = $_;
	}
	$projectData;
}

sub get_seq_data {
	my $file = shift or die "No file name passed to sub";

	my $in = Bio::SeqIO->new(-file => $file, -format => 'genbank');

	my ($top,$alphabet,$sequence,@allfeats,@total_intergenes);
	my ($total_intergenic_len,$contignum) = 0;

	while ( my $seqobj = $in->next_seq ) {
		my @feats = $seqobj->get_SeqFeatures;
		push @allfeats,@feats;
		$sequence .= $seqobj->seq;
		$top = get_topology($seqobj);
		$alphabet = $seqobj->alphabet;
		$contignum++;
		my @arr = get_intergenic_lens(\@feats);
		push @total_intergenes,@arr;
	}

   my $num_intergenic_regions = @total_intergenes;
	( $total_intergenic_len += $_ ) for @total_intergenes;
	my $avg_intergenic_len = sprintf("%d", $total_intergenic_len/$num_intergenic_regions);

	my $pct_gc = get_gc($sequence);
	my $pct_coding = get_coding($sequence,\@allfeats);
	my $genes = grep { $_->primary_tag eq 'gene' } @allfeats;
	my $ncrnas = grep { $_->primary_tag =~ /RNA/ } @allfeats;
	my $avg_gene_len = get_avg_gene_length(\@allfeats);
	my $overlaps = get_overlaps(\@allfeats);
	my $length = length($sequence);

	($pct_gc,$pct_coding,$contignum,$genes,$ncrnas,$avg_gene_len,
	 $avg_intergenic_len,$overlaps,$length,$top,$alphabet);
}

sub get_topology {
	my $seq = shift;
	if ($seq->is_circular) {
		return 'circular';
	} else {
		return 'linear';
	}
}

sub get_gc {
	my $seq = shift;
	my @gcs = $seq =~ /G|C|g|c/g;
	my $pct = (scalar @gcs)/length($seq);
	($pct) = $pct =~ /\d\.(\d\d)\d*/;
	return ($pct);
}

sub get_coding {
	my $seq = shift;
	my $features = shift;
	my @genes = grep { $_->primary_tag eq 'gene' } @$features;
	my $total_cdslength;
	for (@genes) {
		$total_cdslength += $_->length;
	}
	my $pct = $total_cdslength/(length($seq));
	($pct) = $pct =~ /\d\.(\d\d)\d*/;
	return ($pct);
}

sub get_hypothetical {
	my $features = shift;
	my $hypo;
	for (@$features) {
		my ($product) = $_->get_tag_values('product') if ($_->has_tag('product'));
		$hypo++ if ($product && $product =~ /hypothetical/);
	}
	return $hypo;
}

sub get_avg_gene_length {
	my $features = shift;
	my @genes = grep { $_->primary_tag eq 'gene' } @$features;
	my $total_length = 0;
	$total_length += $_->length for (@genes);
	my $avg_gene_length = $total_length / scalar @genes;
	($avg_gene_length) = $avg_gene_length =~ /(\d+)\.*\d*/;
	return $avg_gene_length;
}

sub get_overlaps {
	my $features = shift;
	my @genes = grep { $_->primary_tag eq 'gene' } @$features;
	my ($prev_gene, $num_overlaps);
	$num_overlaps = 0;

	for my $gene (@genes) {
		if ( $prev_gene && $gene->location->overlaps($prev_gene->location) ) {
			$num_overlaps++;
		}
		$prev_gene = $gene;
	}
	return ($num_overlaps);
}

sub get_intergenic_lens {
	my $features = shift;
	my @genes = grep { $_->primary_tag eq 'gene' } @$features;
	my ($prev_gene, @intergenic_length_array);

	for my $gene (@genes) {
		if ( $prev_gene && ! $gene->location->overlaps($prev_gene->location) ) {
			my $intergenic_length = $gene->start - $prev_gene->end;
			push @intergenic_length_array, $intergenic_length;
		}
		$prev_gene = $gene;
	}

	@intergenic_length_array;
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


