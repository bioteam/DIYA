#!/usr/bin/perl -w
# $Id: extract-tar.pl 192 2009-01-16 14:37:11Z briano $

use strict;

my $dest = '/massive/data/genome/submissions';

my $tar = shift or die "No tar file given";

my ($id) = $tar =~ /(\w+)\.tar$/ or die "Could not get id from $tar";

system "cp $id.discrp $dest/$id" if  ( -e "$id.discrp" );
system "cp $id.qual $dest/$id" if  ( -e "$id.qual" );

system "tar xvf $tar";

chdir "$id-gbsubmit.out.d" or die "Cannot chdir to $id-gbsubmit.out.d";

unless ( -e "$dest/$id" ) {
	system "mkdir $dest/$id";
}

for my $suffix ( qw(val gbf qvl sqn tbl fsa agp) ) {
	if ( -e "$id.$suffix" ) {
		system "cp $id.$suffix $dest/$id";
		print "Copied $id.$suffix to $dest/$id\n";
	}
}

system "cp discrp $dest/$id" if  ( -e "discrp" );

__END__

drwxrwxr-x bosborne/bosborne 0 2008-12-04 10:12 yberc0001-gbsubmit.out.d/
-rw-rw-r-- bosborne/bosborne 43298 2008-12-04 10:11 yberc0001-gbsubmit.out.d/discrp
-rw-rw-r-- bosborne/bosborne 4446787 2008-12-04 10:11 yberc0001-gbsubmit.out.d/yberc0001.fsa
-rw-rw-r-- bosborne/bosborne 9831700 2008-12-04 10:11 yberc0001-gbsubmit.out.d/yberc0001.gbf
-rw-rw-r-- bosborne/bosborne 13050947 2008-12-04 10:11 yberc0001-gbsubmit.out.d/yberc0001.qvl
-rw-rw-r-- bosborne/bosborne 26280665 2008-12-04 10:11 yberc0001-gbsubmit.out.d/yberc0001.sqn
-rw-rw-r-- bosborne/bosborne  1090600 2008-12-04 10:11 yberc0001-gbsubmit.out.d/yberc0001.tbl
-rw-rw-r-- bosborne/bosborne     3608 2008-12-04 10:11 yberc0001-gbsubmit.out.d/yberc0001.val
