#!/usr/bin/perl -w

use strict;

use diya::MARC::GenbankConvertUtil;

my $file = shift or die "No file given";

my $gb = diya::MARC::GenbankConvertUtil->new(-file => $file);

$gb->list_primary_tags;
