#!/usr/bin/perl -w
# $Id: rename-and-formatdb.pl 164 2008-08-06 15:39:20Z briano $

=head1 NAME

rename-and-formatdb.pl

=head1 DESCRIPTION

Creates a new file using an existing file name and a strain name 
and runs formatdb on the new file.

=cut

use strict;
use File::Basename qw(fileparse);

my ($oldfullname,$strain) = @ARGV;

die "$0: need both file name and strain name" unless ( $strain && $oldfullname );
die "$0: file $oldfullname" unless ( -e $oldfullname );

my ($filename,$dir,$suffix) = fileparse($oldfullname);

my $newfullname = $dir . $strain . ".fa";

print "Copying $oldfullname to $newfullname\n";
system "cp $oldfullname $newfullname";

print "Running formatdb on $newfullname\n";
system "formatdb -p F -i $newfullname";

if ( -e "formatdb.log" ) {
	print "Removing formatdb.log\n";
	unlink "formatdb.log";
}
