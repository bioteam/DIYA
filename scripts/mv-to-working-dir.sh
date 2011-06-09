#!/bin/bash
# $Id: mv-to-working-dir.sh 233 2009-04-17 17:37:12Z briano $

ARGS=2

if [ $# -ne $ARGS ]  # Correct number of arguments passed to script?
then
  echo "Usage: $0 File Directory"
  exit
fi

cp $1 $2
