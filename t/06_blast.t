#-*-Perl-*-
# $Id: 06_blast.t 193 2008-08-12 16:20:24Z briano $
 
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
my $fa = "little-buchnera.fa";
my $out = "blastp.out";
my $in = "t/data/contig1.fa";

# test command-line options - the blastall command in t/data/blast.conf is
#     <command>-p blastn -d DB -i INPUT -o OUTPUT</command>

SKIP: {

	skip("Skipping test, /usr/local/bin/blastall not found", 2) 
	  unless (-e "/usr/local/bin/blastall");

	system "t/data/diya-blast.pl --set DB=$tmp/$fa --set OUTPUT=$tmp/$out --set INPUT=$in";

	ok( -e "$tmp/$fa.nin", "*nin file created by formatdb" )
	  ||  diag("*nin file not created by formatdb");

	ok( -e "$tmp/$out", "blastn output file created by blastall" )
	  ||  diag("blastn output file not created by blastall");

}

END {
	unlink "formatdb.log";
	system "rm -fr $tmp";
}

1;

__END__
