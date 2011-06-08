#-*-Perl-*-
# $Id: 09_blastpCDS.t 271 2008-12-09 14:05:39Z briano $
 
use strict;
use vars qw($NTESTS);
use lib "./lib";

BEGIN: {
    $NTESTS = 4;
    
    eval { require Test::More;
		     require Test::Exception; };
    if ( $@ ) {
		 use lib 't/lib';
    }
    use Test::More;
	 use Test::Exception;
    
    plan tests => $NTESTS;    
}

my $tmp = "t/tmp/blastCDS";
my $fa = "mini-yberc0001.gbk";
my $t = "t/data";
my $db  = "t/data/example-db.fa";
my $conf = "t/data/blastpCDS.conf";

use_ok('Bio::SeqIO');

SKIP: {

	skip("Skipping test, /usr/local/bin/blastall not found", 3) 
	  unless (-e "/usr/local/bin/blastall");

	skip("Skipping test, /usr/local/bin/formatdb not found", 3) 
	  unless (-e "/usr/local/bin/formatdb");

	# Test a pipeline that runs blasts and annotates a Genbank file
	system "t/data/diya.pl --conf $conf --outputdir $tmp --set MYDB=$db  $t/$fa";

	ok( -e "$db.nin", "*nin file created by formatdb" )
	  ||  diag("*nin file not created by formatdb");

	my @gbks = <$tmp/*out.gbk>;

	is( scalar @gbks, 1, "Genbank file created by blastpCDS" )
	  ||  diag("Genbank file not created by blastpCDS");
	
	my $in = Bio::SeqIO->new(-file => $gbks[0], -format => 'genbank' );
	my $gbk = $in->next_seq;
	my @feats = $gbk->get_SeqFeatures;

	my @score = $feats[1]->get_tag_values('score');

	is( $score[0], 1041, "Score of first CDS is correct" )
	  ||  diag("Score of first CDS is not correct");

}

END {
	unlink "formatdb.log";
	unlink "$db.nin";
	unlink "$db.nsq";
	unlink "$db.nhr";
	system "rm -fr $tmp";
}

1;

__END__

     CDS             complement(267..791)
                     /locus_tag="yberc0001_10"
                     /score="1041"
                     /product="CDS from mini-yberc0001.gbk"
