#!/usr/bin/perl -w
#-------------------------------------------------------------------------------
# Copyright 2008
#
# This file is part of DIYA.
#
#    DIYA is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    DIYA is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with diya.  If not, see <http://www.gnu.org/licenses/>.
#-------------------------------------------------------------------------------

=head1 NAME

Do It Yourself Annotator - diya.pl

=head1 DESCRIPTION

This is a template for a script to run a diya pipeline.

=head1 COMMAND-LINE OPTIONS

=over 5

=item B<--verbose>

Set the verbosity level, 0 or 1.

 diya.pl --verbose 1

=item B<--inputfile>

 diya.pl --inputfile B1000.seq

=item B<--mode> [serial|parallel]

Run the batch in serial mode, or parallel mode if SGE is available.

 diya.pl --mode serial

=item B<--set> tag=value

Manually set a configuration option. To see the available options do:

 diya.pl --set

To change or set an option do something like:

 diya.pl --set 

=item B<--save> 

If the configuration has changed then you can save this new
configuration to a file. You have to provide a file name. For example:

 diya.pl --save new-diya.conf --mode parallel

=back	

=cut

use strict;
use lib "./lib";
use diya;

my $diya = diya->new( -use_conf  => "t/data/blast.conf",
		      -verbose   => 1 );

$diya->inputfile("t/data/little-buchnera.fa");
$diya->outputdir("t/tmp/Blast");
$diya->read_conf;
$diya->run;
