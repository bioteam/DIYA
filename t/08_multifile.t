#-*-Perl-*-
# $Id: 08_multifile.t 195 2008-08-12 16:24:24Z briano $
 
use strict;
use vars qw($NTESTS);
use lib "./lib";

BEGIN: {
    $NTESTS = 2;
    
    eval { require Test::More;
		     require Test::Exception; };
    if ( $@ ) {
        use lib 't/lib';
    }
    use Test::More;
	 use Test::Exception;

    plan tests => $NTESTS;    
}

my $tmp = "t/tmp/Blast";
my $fa1 = "t/data/little-buchnera.fa";
my $fa2 = "t/data/contig1.fa";
my $db  = "t/data/little-buchnera.fa";
my $conf = "t/data/blast2.conf";

# Test that multiple input files are processed correctly,
# fed one after the other to the pipeline

SKIP: {

	skip("Skipping test, /usr/local/bin/blastall not found", 2) 
	  unless (-e "/usr/local/bin/blastall");

	skip("Skipping test, /usr/local/bin/formatdb not found", 2) 
	  unless (-e "/usr/local/bin/formatdb");

	system "t/data/diya-blast2.pl --conf $conf --set MYDB=$db $fa1 $fa2";

	ok( -e "$db.nin", "*nin file created by formatdb" )
	  ||  diag("*nin file not created by formatdb");

	my @blastouts = <$tmp/*-blastall.out>;
	print "blastall output files are @blastouts\n";

	is( $#blastouts, 1, "2 blastn output files created by blastall" )
	  ||  diag("2 blastn output files not created by blastall");

}

END {
	unlink "formatdb.log";
	system "rm -fr $tmp";
	unlink "$db.n*";
}

1;

__END__
