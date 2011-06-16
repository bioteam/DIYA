#!/bin/bash

ARGS=2
DEST=/massive/data/genome/submissions
DIR=$1
ID=$2

if [ $# -ne $ARGS ]  # Correct number of arguments passed to script?
then
  echo "Usage: $0 Directory Id"
  exit
fi

if [ ! -d $DEST/$ID ]
then
   mkdir $DEST/$ID
	echo "Making directory $DEST/$ID"
fi

if [ -d $DIR ]
then

   for f in $ID.*
	do
   	mv $f $DIR
	done

   cd $DIR

	for f in $( ls *gbk *out $ID* 454AllContigs.fna )
	do
   	cp $f $DEST/$ID
		echo "Copying $f to $DEST/$ID"
	done

else
	echo "No directory $DIR found"
fi

cd ..

for f in $( ls *err *out )
do
   rm $f
done

if [ -d $ID-gbsubmit.out.d ]
then
	cd $ID-gbsubmit.out.d

   for f in $( ls $ID* discrp )
   do
      cp $f $DEST/$ID
		echo "Copying $f to $DEST/$ID"
   done

else
	echo "No directory $ID-gbsubmit.out.d found"
fi

echo "Done with copy and cleanup for $ID and $DIR"