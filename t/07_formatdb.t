#-*-Perl-*-
# $Id: 07_formatdb.t 270 2008-12-09 14:02:48Z briano $
 
use strict;
use vars qw($NTESTS);
use lib "./lib";

BEGIN: {
    $NTESTS = 7;
    
    eval { require Test::More;
		     require Test::Exception; };
    if ( $@ ) {
        use lib 't/lib';
    }
    use Test::More;
	 use Test::Exception;

    plan tests => $NTESTS;    
}

my $testconf = "t/data/blast.conf";
my $tmpdir = "t/tmp/blast";
my $script = 'formatdb';
my $module = 'blastall';

use_ok('diya');

my $blast = diya->new();
isa_ok ( $blast, "diya" );
$blast->verbose(1);
$blast->read_conf($testconf);
$blast->outputdir($tmpdir);
$blast->inputfile("t/data/little-buchnera.fa");

my $type = $blast->_get_type($script);
is($type, 'script', "$script is of type script");

$type = $blast->_get_type($module);
is($type, 'parser', "$module is of type parser");
#
# check <inputfrom>
#
my $inputfrom = $blast->_inputfrom("blastall");
is($inputfrom,'formatdb',"Can get inputfrom() of parser from $testconf" || 
  "Can not get inputfrom() from $testconf");

$blast->order('formatdb');

SKIP: {

	skip("Skipping test, /usr/local/bin/formatdb not found", 2) 
	  unless (-e "/usr/local/bin/formatdb");

	$blast->run;

	ok(-e "$tmpdir/little-buchnera.fa.nin", "*nin file created by formatdb" )
	  ||  diag("*nin file not created by formatdb");

	$blast->cleanup;
	ok(! -d $tmpdir, "Directory $tmpdir removed by cleanup()" ) 
	  || diag("Directory $tmpdir not removed by cleanup()");

}

1;

END {
	unlink "formatdb.log";
}

__END__
