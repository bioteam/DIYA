#!/usr/bin/perl
# Make individual *cm files from *aln files

use strict;

my @alns = <aln/*.aln>;

for my $aln (@alns) {

    my $line = `grep 'GF AC' $aln`;
    my ($id) = $line =~ /(RF\d+)/;

    my $cmd = "/usr/local/share/infernal-1.0.1/bin/cmbuild $id.cm $aln";
    print "Command is $cmd\n";
    `$cmd`;
    #$cmd = "/usr/local/share/infernal-1.0.1/bin/cmcalibrate $id.cm";
    #print "Command is $cmd\n";
    #`$cmd`;
    print "Done with $1.cm\n";

}

__END__

