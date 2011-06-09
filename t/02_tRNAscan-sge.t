#-*-Perl-*-
 
use strict;
use vars qw($NTESTS);
use lib "./lib";

BEGIN: {
    $NTESTS = 6;
    
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
my $tmpdir   = "t/tmp/trnaoutputsge";
my $module   = 'MARC::tRNAscanSE';

use_ok('diya');

my $scanner = diya->new(-mode => "sge");
isa_ok ( $scanner, "diya" );
$scanner->verbose(1);
$scanner->read_conf($testconf);
$scanner->inputfile($gbk_file);
$scanner->outputdir($tmpdir);
is ($scanner->mode(), "sge", "submission mode is not sge.");

my $scan = "/usr/local/bin/tRNAscan-SE";

SKIP: {

	skip("Can't run tests, $scan not found", 2) 
	  unless (-e "$scan");

	skip('Skipping SGE tests, $SGE_ROOT variable not set', 2)
	  unless ( defined $ENV{SGE_ROOT} );

	$scanner->run;

	ok(defined($scanner->_lastsgeid()), "job is submitted to SGE" )
	  ||  diag("SGE submisson fails");
  
	ok(-e "$tmpdir/BA000003-nt.fasta", "Fasta file created from Genbank file" )
	  ||  diag("Fasta file not created from Genbank file");

}

$scanner->cleanup;
ok(! -d $tmpdir, "Directory $tmpdir removed by cleanup()" ) 
  || diag("Directory $tmpdir not removed by cleanup()");


1;

__END__
