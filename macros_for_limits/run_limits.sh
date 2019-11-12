#!/bin/bash

LISTDIR="list_datacards.txt"

i=0
j=0
#Cleaning folders

for run in `cat $LISTDIR`

do

    i=$((i+1))
    echo "cat $LISTDIR_outputname/list_outnames_ossf/${i}"
    echo "****> ${run} "

    JOBNAME_out="$1/"

    echo "Analizzando run ${run}_datacard.root"

    echo "Lanciando job #${i}"
    echo "jobout ${JOBNAME_out}_datacard.root"
    echo "combine -M AsymptoticLimits --name ${run} ${run}_datacard.root"

    combine -M AsymptoticLimits --name ${run} ${run}_datacard.root

done
