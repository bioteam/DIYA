#!/bin/bash

formatdb -p $2 -V T -i $1

if [[ -e formatdb.log ]] 
then
	rm formatdb.log
fi

