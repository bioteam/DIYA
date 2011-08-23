#!/bin/bash

ARGS=2

if [ $# -ne $ARGS ]  # Correct number of arguments passed to script?
then
  echo "Usage: $0 File Directory"
  exit
fi

ID=$1
DIR=$2

mv -f $ID.* phage_phinder_${ID}.log $DIR
rm -fr *err *out GLOCAL* FRAG* combined* error.log ${ID}_dir phage_finder*
