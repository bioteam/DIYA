#!/bin/bash

ARGS=2
INPUT=$2
ID=$1

if [ $# -ne $ARGS ]  # Correct number of arguments passed to script?
then
  echo "Usage: $0 Id Input"
  exit
fi

rm $ID.gbk

cp $INPUT.gbk $ID.gbk
