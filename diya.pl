#!/arch/bin/perl
#-------------------------------------------------------------------------------
# Copyright 2008
#
# This file is part of DIYA.
#
# DIYA is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# DIYA is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with the diya package.  If not, see <http://www.gnu.org/licenses/>.
#--------------------------------------------------------------------------

=head1 NAME

Do It Yourself Annotator - diya.pl

=head1 SYNOPSIS

diya.pl [options] [inputfiles]

=head1 DESCRIPTION

This is a template for a script to run a diya pipeline. For more
information see the diya documentation in docs/diya.html or do:

 >perldoc diya

=cut

use strict;
use warnings;
use diya;

my $diya = diya->new( -verbose => 1);

$diya->read_conf;

$diya->run;
