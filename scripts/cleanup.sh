#!/bin/bash
# $Id: cleanup.sh 256 2009-05-14 01:30:41Z briano $

ARGS=2

if [ $# -ne $ARGS ]  # Correct number of arguments passed to script?
then
  echo "Usage: $0 File Directory"
  exit
fi

DIR=$1
ID=$2

mv $ID.* $DIR
rm *err *out
