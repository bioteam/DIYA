#!/bin/bash
# $Id: formatdb.sh 234 2009-04-17 17:58:35Z briano $

ARGS=2
IN=$1
BOOL=$2
FORMATDB=/usr/local/bin/formatdb

if [ $# -ne $ARGS ]  # Correct number of arguments passed to script?
then
  echo "Usage: $0 Input T/F"
  exit
fi

if [ -x $FORMATDB ]
then
	$FORMATDB -p $BOOL -V T -i $IN
else
	echo "$FORMATDB not found"
	exit
fi

if [ -e formatdb.log ]
then
	rm formatdb.log
fi