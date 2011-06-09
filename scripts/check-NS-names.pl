#!/usr/bin/perl
# $Id: check-NS-names.pl 217 2009-02-27 15:36:48Z briano $

=head1 NAME

check-NS-names.pl

=head1 DESCRIPTION

Gets the correct strain names from NCBI, checks them against the 'gbtaxname'
field for all NS pages.

=cut

use strict;
use Perlwikipedia;
use Bio::DB::EUtilities;
use XML::Simple;

my $wiki = { 'host' => 'loki.bdrd',
				 'dir'  => 'wiki',
				 'user' => 'Bioteam',
				 'pass' => 'bioteampw' };

my $bot = Perlwikipedia->new;
$bot->set_wiki($wiki->{host}, $wiki->{dir}); 
$bot->login($wiki->{user}, $wiki->{pass});

my @strains = $bot->get_pages_in_category('Category:Strain');

for my $strain ( @strains ) {

	my $text = $bot->get_text($strain);
	next if ( $text == '2' );

	my ($taxid) = $text =~ /\s*\|\s*taxid\s*=\s*(\d+)/;
	next unless $taxid;

	my $new_name = get_name($taxid);
	next unless $new_name;

	my ($old_name) = $text =~ /\s*\|\s*gbtaxname\s*=\s*([^\n]*)\s*/;

	if ( $new_name ne $old_name ) {
		print "Strain:$strain\tOld:$old_name\tNew:$new_name\n";
	}
}

print "Done\n";

sub get_name {
	my $id = shift;

   my $factory = Bio::DB::EUtilities->new(-eutil => 'efetch',
                                          -db    => 'taxonomy',
                                          -id    => $id );

   my $res = $factory->get_Response->content;

	my $data = XMLin($res);

	$data->{Taxon}->{ScientificName};
}

__END__


