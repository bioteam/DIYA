#-*-Perl-*-

use strict;
use vars qw($NTESTS);
use lib "./lib";

BEGIN: {
    $NTESTS = 8;

    eval {
        require Test::More;
        require Test::Exception;
    };
    if ($@) {
        use lib 't/lib';
    }
    use Test::More;
    use Test::Exception;

    plan tests => $NTESTS;
}

my $gbk      = 't/data/2012_03_11_14_27_03-MARC::cmsearch.out.gbk';
my $tbl2asn  = '/usr/local/bin/tbl2asn';
my $script   = 'scripts/gbconvert.pl';
my $template = 'examples/template';
my $gdir     = 'test0001-gbsubmit.out.d';

use_ok('diya');

SKIP: {

    skip( "Skipping test, file $gbk not found", 7 )
      unless ( -e $gbk );
    skip( "Skipping test, script $script not found", 7 )
      unless ( -x $script );
    skip( "Skipping test, file $template not found", 7 )
      unless ( -e $template );
    skip( "Skipping test, application $tbl2asn not found", 7 )
      unless ( -x $tbl2asn );

	`$script -template $template -e $tbl2asn -a test $gbk`;

    ok( -s "$gdir/test0001.fsa", "*fsa file created by $script" )
      || diag("*fsa file not created by $script");
    ok( -s "$gdir/test0001.gbf", "*gbf file created by $script" )
      || diag("*gbf file not created by $script");
    ok( -s "$gdir/test0001.sqn", "*sqn file created by $script" )
      || diag("*sqn file not created by $script");
    ok( -s "$gdir/test0001.tbl", "*tbl file created by $script" )
      || diag("*tbl file not created by $script");
    ok( -s "$gdir/discrp", "discrp file created by $script" )
      || diag("discrp file not created by $script");
    ok( -z "$gdir/test0001.val", "empty *val file created by $script" )
      || diag("*val file is not empty");
    ok( -z "$gdir/errorsummary.val", "empty errorsummary.val file created by $script" )
      || diag("errorsummary.val file is not empty");
}

1;

END {
    `rm -fr $gdir`;
}

__END__
