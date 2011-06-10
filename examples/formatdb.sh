#!/bin/csh
# $Id: formatdb.sh 172 2008-08-06 20:24:47Z briano $

/usr/local/share/apps/ncbi/bin/formatdb -p $argv[2] -V T -i $argv[1]

rm formatdb.log
