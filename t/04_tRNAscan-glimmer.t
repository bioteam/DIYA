#-*-Perl-*-
# $Id: 04_tRNAscan-glimmer.t 192 2008-08-12 16:14:26Z briano $
 
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

my $gbk_file = "t/data/BA000003-nt.gbk";
my $testconf = "t/data/tRNAscan.conf";
my $tmpdir = "t/tmp/trna-glimmer-out";

use_ok('diya');

my $scanner = diya->new();
$scanner->verbose(1);
$scanner->read_conf($testconf);
$scanner->inputfile($gbk_file);
$scanner->outputdir($tmpdir);

my @order = qw( tRNAscanSE glimmer3 );
$scanner->order(@order);

SKIP: {

	for my $step (@order) {
		my $path = $scanner->_home($step) . $scanner->_executable($step);
		skip("Skipping run() step, $path not found",2) unless (-e $path);
	}

	$scanner->run;

	# check for glimmer3 file and Genbank file

	my $out = $scanner->_outputfile('glimmer3') . ".detail";
	my $gbk = $scanner->_outputfile('glimmer3') . ".gbk";

	ok(-e $out, "$order[1] output file created" ) 
	  || diag("$order[1] output file not created");

	ok(-e $gbk, "$gbk output file for $order[1] created" ) 
	  || diag("Genbank output file for $order[1] not created");

}

$scanner->cleanup;
ok(! -d $tmpdir, "Directory $tmpdir removed by cleanup()" ) 
  || diag("Directory $tmpdir not removed by cleanup()");


1;

__END__
