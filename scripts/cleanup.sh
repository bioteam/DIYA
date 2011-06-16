#!/bin/bash

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
