#!/usr/bin/perl
# $Id: run-diya-on-project.pl 268 2009-07-14 23:40:37Z briano $

=head1 NAME

compare-2-genomes.pl

=head1 DESCRIPTION

Script to compare 2 genomes annotated by Diya. Example:

 ./compare-2-genomes.pl -s banth007 bthur004

Run this script in a directory containing the Diya results directories.

=cut

use strict;
use Bio::SeqIO;
use Bio::Tools::Run::StandAloneBlast;
use Getopt::Long;
use File::Find;

my (@strains,%strains);
my $threshold = '100';
my $verbose = 0;

GetOptions( 's=s{2,}' => \@strains, "t=i" => \$threshold, "verbose" => \$verbose );

die "Must supply 2 strain names with \'-s\'" if ( $#strains != 1 );

for my $strain ( @strains ) {

	die "Invalid strain name" if ( $strain !~ /^[a-zA-Z]{5}\d{4}$/ );

	my $gbk = findFile($strain);

	my $in = Bio::SeqIO->new( -file => $gbk, -format => 'genbank' );
	while ( my $seq = $in->next_seq ) {

		for my $tag ( qw( CDS gene tRNA rRNA ) ) {

			my @feats = grep { $_->primary_tag eq $tag } $seq->get_SeqFeatures;
			push @{$strains{$strain}{$tag}}, @feats;

		}
	}
}

# Create fasta files
for my $strain ( keys %strains ) {
	for my $tag ( qw( CDS gene tRNA rRNA ) ) {

		`rm $strain-$tag.fa*` if ( -e "$strain-$tag.fa" );
		my $out = Bio::SeqIO->new(-file => ">>$strain-$tag.fa",-format => 'fasta');

		for my $feat ( @{$strains{$strain}{$tag}}  ) {
			if ( $tag eq 'CDS' ) {
				my @peps = $feat->get_tag_values('translation') if ( $feat->has_tag('translation') );
				my @ids = $feat->get_tag_values('locus_tag') if ( $feat->has_tag('locus_tag') );
				my $pep = Bio::Seq->new(-seq => $peps[0], -display_id => $ids[0] );
				$out->write_seq($pep);
			} else {
				my @ids = $feat->get_tag_values('locus_tag') if ( $feat->has_tag('locus_tag') );
				my $seq = $feat->seq;
				$seq->display_id($ids[0]);
				$out->write_seq($seq);
			}
		}
		$tag eq 'CDS' ? `formatdb -p T -i $strain-$tag.fa` : `formatdb -p F -i $strain-$tag.fa`;
	}
}

for my $strain ( keys %strains ) {

 OTHER:
	for my $otherStrain ( keys %strains ) {

		next OTHER if ( $strain eq $otherStrain );

		for my $tag ( qw( CDS gene tRNA rRNA ) ) {

			my $program = 'blastn';
			$program = 'blastp' if ( $tag eq 'CDS' );

			print "$strain has " . scalar @{$strains{$strain}{$tag}} . " $tag features\n";

		  	my $out = Bio::SeqIO->new(-file => ">>$tag-unique-$threshold-$strain.fa", -format => 'fasta' );
			my $count;
			my $all;

			# Create BLAST database
			my @params = (program  => $program,
							  database => "$otherStrain-$tag.fa");
			my $factory = Bio::Tools::Run::StandAloneBlast->new(@params);

		 FEAT:
			for my $feat ( @{$strains{$strain}{$tag}} ) {
				my $seq = $feat->seq;
				my @ids = $feat->get_tag_values('locus_tag') if ( $feat->has_tag('locus_tag') );
				$seq->display_id($ids[0]);

				if ( $tag eq 'CDS' ) {
					my @peps = $feat->get_tag_values('translation') if ( $feat->has_tag('translation') );
					$seq->seq($peps[0]);
				}

				my $len = $seq->length;
				my $report = $factory->blastall($seq);
				if ( isUnique($report,$threshold,$len) ) {
					$out->write_seq($seq);
					$count++;
					print "Total unique $tag ($strain): $count\n" if $verbose;
				}
				$all++;
				print "Total $tag ($strain): $all\n" if $verbose;
			}
			print "$strain has $count unique $tag features (threshold $threshold)\n";
		}
	}
}

sub isUnique {
	my ($report,$threshold,$len) = @_;

	while ( my $result = $report->next_result ) {
		while ( my $hit = $result->next_hit ) {
			while ( my $hsp = $hit->next_hsp ) {
				return 0 if ( $hsp->length('hit') == $len && $hsp->percent_identity >= $threshold );
			}
      }
    }
  1;
}

sub findFile {
	my $name = shift;
	my @files;

	find( sub { push @files, $File::Find::name if -f && /^$name.gbf$/ }, "$name-gbsubmit.out.d" );

	die "More than one $name file found in $name-gbsubmit.out.d" if ( $#files > 0 );
	die "No $name file found in $name-gbsubmit.out.d" if ( ! @files );

	$files[0];
}

__END__


