#!/usr/bin/perl
# Make individual *aln files from Rfam.full

use strict;

open MYIN,"Rfam.full";
my $num = 0;

while (my $line = <MYIN>) {

    if ( $line =~ /^# STOCKHOLM/ ) {
	 $num++;
	 my $id = 'RF' . sprintf("%05d", $num);
	 open MYALN,">>aln/$id.aln";
    }

    print MYALN $line;

}

__END__

    `/usr/local/share/infernal-1.0.1/bin/cmbuild $id.cm aln/$id.aln`;
    `/usr/local/share/infernal-1.0.1/bin/cmcalibrate $id.cm`;
    print "Done with $id.cm\n";
