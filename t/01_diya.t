#-*-Perl-*-
 
use strict;
use vars qw( $NTESTS );
use lib "./lib";

BEGIN: {
	$NTESTS = 48;

    eval { require Test::More;
		     require Test::Exception; };
    if ( $@ ) {
        use lib 't/lib';
    }
    use Test::More;
	 use Test::Exception;

	plan tests => $NTESTS;
}

use_ok("diya");
use_ok("diya::MARC::CRT");
use_ok("diya::MARC::rnammer");
use_ok("diya::MARC::rpsblastCDS");
use_ok("diya::MARC::GenbankConvertUtil");
use_ok("diya::MARC::phobos");
use_ok("diya::MARC::blastpCDS");
use_ok("diya::MARC::glimmer3");
use_ok("diya::MARC::blastall");
use_ok("diya::MARC::tRNAscanSE");
use_ok("diya::MARC::cmsearch");

my $diya = diya->new();
isa_ok($diya,'diya');
#
# Set or get various values
#
$diya->verbose(1);
is($diya->verbose, 1, "Can set verbose to 1" || "Cannot set verbose to 1");
$diya->verbose(0);
is($diya->verbose,'0',"Can set verbose to 0");
$diya->verbose(1);

$diya->_use_conf("imaginary");
is($diya->_use_conf,"imaginary","Can set _use_conf" || "Cannot set _use_conf");

$diya->inputfile("imaginary");
is($diya->_next_inputfile, "imaginary","Can get input file" || "Cannot get input file");

if (defined $ENV{DIYAHOME}) {
	my $confPath = $ENV{DIYAHOME} . "/diya.conf";
	ok(-e $confPath, '_diyahome works when $DIYAHOME is set') ||
	  diag("Could not find diya.conf in $ENV{DIYAHOME}");
} else {
	my $pwd = `pwd`;
	chomp $pwd;
	is($diya->_diyahome, "$pwd/", '_diyahome works when $DIYAHOME is not set');
}
#
# set values using new()
# 
my $tmpdir = "t/tmp/test/";

my $p = diya->new(-verbose   => 1,
						-use_conf  => "conf",
						-outputdir => "$tmpdir",
					  );
is($p->_use_conf, "conf", "Can set use_conf with new()") || diag("Cannot set use_conf");
is($p->verbose, 1, "Can set verbose with new()") || diag("Cannot set verbose");
is($p->outputdir, $tmpdir, "Can set outputdir with new()") || diag("Can not set outputdir with new");
#
# Read diya.conf and check
#
my $testconf = "t/data/test.conf";
$diya->_use_conf($testconf);
$diya->read_conf;
ok(defined $diya->_conf, "Can read $testconf") or diag("Cannot read $testconf");
my $conf = $diya->_conf;
check_keys($conf);

sub check_keys {
	my $conf = shift;
	foreach my $name ( qw( run script parser order ) ) {
		my $res = $conf->first_child($name);
		ok( defined $res, "Found <$name> in $testconf") or 
		  diag("Did not find <$name> in $testconf");
	}
}
# TODO
#
# Round-trip a conf file using write_conf()
#
#my $tmpconf = "t/tmp/testconf";
# $diya->write_conf($tmpconf);
# ok(-e $tmpconf, "Can write $tmpconf") || diag("$tmpconf not written");
# $diya->read_conf($tmpconf);
# my $newconf = $diya->_conf;
# ok(defined $newconf, "Can re-read conf") or diag("Cannot re-read conf file");
# is_deeply($conf,$newconf,"Can round-trip a conf file") or 
#   diag("Cannot round-trip a conf file");
#
# make a conf file when no name is supplied
#
# my $confname = $diya->write_conf;
# ok(-e $confname, "Can write conf file with timestamp") || 
#   diag("Cannot write conf file with timestamp");
#
# get and set order() 
#
$diya->_use_conf($testconf);
$diya->read_conf;
my @modules = $diya->order;
my $mod_str = join " ",@modules;
is($mod_str,"MARC::tRNAscanSE MARC::glimmer3 MARC::blastall","Can get with order()") ||
  diag("Can not get with order()");

$diya->order("MARC::blastall");
@modules = $diya->order();
is($modules[0],"MARC::blastall","Can set single step with order()") ||
  diag("Can not set single step with order()");

my @new_modules = qw(MARC::glimmer3  MARC::blastall); 
$diya->order( @new_modules );
@modules = $diya->order();
is_deeply(\@modules, \@new_modules,"Can set order() with array") ||
  diag("Can not set order() with array");

# pass a bad step name
dies_ok{ $diya->order("non-existent") } 'order() dies if passed a bad name';
#
# check _parsers
#
my @ps = qw( MARC::blastall  MARC::glimmer3 MARC::tRNAscanSE );
my @cs = $diya->_parsers;
is_deeply(\@cs, \@ps, "Can get with _parsers" || "Can not get with _parsers");
#
# get executable name for parser
#
$p = $diya->_executable("MARC::blastall");
is($p, "blastall", "Can get executable name with _executable") ||
diag("Can not get executable name with _executable");
#
# get inputformat for parser
#
$p = $diya->_inputformat("MARC::blastall");
is($p, "fasta", "Can get format with _inputformat") ||
diag("Can not get format with _inputformat");
#
# get command for parser
#
$p = $diya->_command("MARC::blastall");
is($p, "-p blastp -d ran.fa -i INPUTFILE > OUTPUTFILE", "Can get command with _command") ||
diag("Can not get command with _command");
#
# get inputfrom for parser
#
$p = $diya->_inputfrom("MARC::blastall");
is($p, '', "Can get inputfrom with _inputfrom") ||
diag("Can not get inputfrom with _inputfrom");
#
# get home for parser
#
$p = $diya->_home("MARC::blastall");
like($p, qr(/usr/local/share/apps/ncbi/bin), "Can get home with _home") ||
diag("Can not get home with _home");
#
# get and set mode
#
$diya->mode("sge");
is($diya->mode, "sge","Can set mode" || "Can not set mode");
# pass a bad mode name
dies_ok{ $diya->mode("non-existent") } 'mode() dies if passed a bad name';
#
# check that output directories are made
#
my $t = diya->new;
$t->verbose(1);
my $timestampdir = $t->_check_outputdir;
ok(-d $timestampdir, "Can create output directory $timestampdir") || 
  diag("Can not create output directory $timestampdir");

my $x = diya->new(-outputdir => $tmpdir);
$x->_check_outputdir;
ok(-d $tmpdir, "Can create output directory $tmpdir") || 
  diag("Can not create output directory $tmpdir");
#
# check that read_conf() does not overwrite values set previously
#
$t = diya->new();
$t->mode("sge");
$t->read_conf($testconf);

is( $t->mode, "sge", "read_conf() does not overwrite existing values" ) || 
  diag("read_conf() overwrites existing values");
#
# check that the type ('script','parser') is retrieved
#
my $type = $t->_get_type("diya-postbatch");
is($type,'script',"Can get type of script from $testconf" || "Can not get type");

$type = $t->_get_type("MARC::blastall");
is($type,'parser',"Can get type of parser from $testconf" || "Can not get type");

$conf = "t/data/blast.conf";
my $z = $diya->new(-use_conf => $conf);
$z->read_conf;

$type = $z->_get_type("formatdb");
is($type,'script',"Can get type of script from $conf" || "Can not get type from $conf");

$type = $z->_get_type("MARC::blastall");
is($type,'parser',"Can get type of parser from $conf" || "Can not get type from $conf");
#
# check <inputfrom>
#
my $inputfrom = $z->_inputfrom("MARC::blastall");
is($inputfrom,'formatdb',"Can get inputfrom() of parser from $conf" || 
  "Can not get inputfrom() from $conf");
$inputfrom = $z->_inputfrom("formatdb");
is($inputfrom, '', "Can get null inputfrom() of parser from $conf" || 
  "Can not get null inputfrom() from $conf");
#
# _check_executable
# 
$conf = "t/data/tRNAscan.conf";
my $q = $diya->new(-use_conf => $conf);
$q->read_conf;
# pass a non-existent application name
dies_ok{ $diya->_check_executable("imaginary") } 
  '_check_executable() dies if passed a bad name';

my $path = 'scripts/g3-from-scratch.csh';

SKIP: {
        skip("Can't run 1 test, $path not found", 1) 
          unless (-e "$path");

		  my $result = $q->_check_executable("glimmer3");
		  is($result, 1,"_check_executable finds executable" || 
			  "_check_executable does not find executable");
	  }


END: {
	#unlink $tmpconf;
	#unlink $confname;
	system "rm -fr $timestampdir";
	system "rm -fr $tmpdir";
}

1;

__END__

