#!/bin/csh
# $Id: cleanup.sh 269 2008-12-04 03:34:05Z briano $

setenv ID $argv[2];
setenv DIR $argv[1];

mv ${ID}.* $DIR;

rm *err;
rm *out;

