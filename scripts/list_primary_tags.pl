#!/usr/bin/perl -w
# $Id: list_primary_tags.pl 235 2009-04-17 19:28:38Z briano $

use strict;

use diya::GenbankConvertUtil;

my $file = shift or die "No file given";

my $gb = diya::GenbankConvertUtil->new(-file => $file);

$gb->list_primary_tags;
