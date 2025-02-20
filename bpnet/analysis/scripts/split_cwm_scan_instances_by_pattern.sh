#!/bin/bash

#Melanie Weilert
#January 2020
#Purpose: Take CWM-scanned list of pattern instances, break into
#separate TSV files by pattern

 #TASKS=$1 #comma separated list of tasks to run through
 #MODISCODIR=$2 #modisco directory up until tasks to parse through.

 # for TASK in $(echo $TASKS | sed "s/,/ /g")
 # do
 #     echo "Separating patterns of $TASK..."

 #     #Define paths
  #   INSTANCES="${TASK}/motif-instances-all-regions.tsv.gz"  
 #    CWMDIR=$(dirname "${INSTANCES}")
 #    NEWDIR=instances_by_pattern
 #    #Create folder to store instances
 #    mkdir $NEWDIR
 #    cd $NEWDIR

     #Keep header
    #gzip -cd ${INSTANCES} | head -1 > header_tmp.txt

    #Unzip | remove header | split by pattern_short column, save with appropriate name
    #gzip -cd $INSTANCES | awk 'NR>1' | awk -F '\t' '{print > (""$26"_motif-instances-all-regions.tsv")}' #separate by pattern

    echo "Adding header and rezipping for $TASK..."

    #add header and zip
    for file in *_motif-instances-all-regions.tsv
     do
      cat header_tmp.txt ${file} | gzip > ${file}.gz
    done

    #remove unzipped files and tmp files
    rm *_motif-instances-all-regions.tsv
    rm header_tmp.txt

#done
