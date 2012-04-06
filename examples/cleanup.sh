#!/bin/csh

setenv ID $argv[2];
setenv DIR $argv[1];

mv ${ID}.* $DIR;

rm *err;
rm *out;

