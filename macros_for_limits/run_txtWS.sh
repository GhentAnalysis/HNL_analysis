#!/bin/bash

LISTDIR="list_datacards.txt"


i=0
j=0
#Cleaning folders


for run in `cat $LISTDIR`

do
    i=$((i+1))

    echo "Analizzando run ${run}"
    echo "Lanciando job #${i}"
    echo "text2workspace.py ${run}"

    ../../../scripts/text2workspace.py  ${run}_datacard.txt

done
