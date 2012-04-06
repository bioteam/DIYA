#!/bin/csh

/usr/local/share/apps/ncbi/bin/formatdb -p $argv[2] -V T -i $argv[1]

rm formatdb.log
