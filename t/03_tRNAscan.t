#-*-Perl-*-
# $Id: 03_tRNAscan.t 191 2008-08-12 16:09:38Z briano $
 
use strict;
use vars qw($NTESTS);
use lib "./lib";

BEGIN: {
    $NTESTS = 11;
    
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
my $tmpdir = "t/tmp/trnaoutput";
my $module = 'tRNAscanSE';

use_ok('diya');

my $scanner = diya->new();
isa_ok ( $scanner, "diya" );
$scanner->verbose(1);
$scanner->read_conf($testconf);
$scanner->inputfile($gbk_file);
$scanner->outputdir($tmpdir);

my $string = "-B -o OUTPUTFILE INPUTFILE";
is($scanner->_command($module), $string, "Can get <command> for $module");

$string = "serial";
is($scanner->mode("serial"), $string, "Can get <mode> for the pipeline");

my $home = "/usr/local/bin/";
is($scanner->_home($module), $home, "Can get <home> for $module");

my $exe = "tRNAscan-SE";
is($scanner->_executable($module), $exe, "Can get <home> for $module");

SKIP: {

	skip("Can't run tests, $home$exe not found", 4) 
	  unless (-e "$home$exe");

	$scanner->run;

	ok(-e "$tmpdir/BA000003-nt.fasta", "Fasta file created from Genbank file" )
	  ||  diag("Fasta file not created from Genbank file");
	
	is($scanner->_next_inputfile, "0", "File queue is empty" 
		|| "File queue is not empty but it should be");

	my $out = $scanner->_outputfile($module);
	ok(-e $out, "$module output file created" ) || diag("$module output file not created");

	my $gbk = $out. ".gbk";
	ok(-e $gbk, "Genbank output file created" ) || diag("Genbank output file not created");

}

$scanner->cleanup;
ok(! -d $tmpdir, "Directory $tmpdir removed by cleanup()" ) 
  || diag("Directory $tmpdir not removed by cleanup()");


1;

__END__
