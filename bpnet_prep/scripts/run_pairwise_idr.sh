#!bin/bash

#Melanie Weilert
#August 2019
#Purpose: Given a string representing a factor, run pairwise (non-repeating) combinations of IDR
#Based on Snakemake file setup.

WDIR=$1 #input working directory
COMBINED_MACS2_FILE=$2 #input factor
cd $WDIR

#Define files
FACTOR=$(echo ${COMBINED_MACS2_FILE} | sed 's/.*mtsc_\(.*\)_nexus.*/\1/')
REPLICATE_MACS2_FILES="macs2/individual/mtsc_${FACTOR}_nexus_*_peaks.narrowPeak"

echo "Running IDR for ${FACTOR}..."
COMBOSRAN=()
#Run pairwise replicates of idr
for i in $REPLICATE_MACS2_FILES; do
  for j in $REPLICATE_MACS2_FILES; do
    if [ "$i" != "$j" ]; then
      #define samples
      REPA=$(echo ${i} | sed 's/.*_nexus_\(.*\)_peaks.*/\1/')
      REPB=$(echo ${j} | sed 's/.*_nexus_\(.*\)_peaks.*/\1/')

      #Record this as a combination ran
      ID1=$REPA\_$REPB
      ID2=$REPB\_$REPA

      #If combos do not contain the ID already
      if [[ ! " ${COMBOSRAN[@]} " =~ " ${ID1} " ]] && [[ ! " ${COMBOSRAN[@]} " =~ " ${ID2} " ]]
      then
        echo "Running ${ID1} and ${ID2}..."
        COMBOSRAN+=($ID1)
        #echo "${COMBOSRAN[@]}"
        OUTPUT="idr/pairwise/${FACTOR}/mtsc_${FACTOR}_${REPA}_vs_${REPB}_nexus_idr.txt"
        #run IDR
        idr --samples $i $j --peak-list $COMBINED_MACS2_FILE \
        --input-file-type narrowPeak --output-file {} \
        --output-file $OUTPUT --idr-threshold 0.05
      fi
    fi
  done
done

# #Once pairwise combinations are done, select pair that has the most returned peaks
# #This is consistent with the ENCODE recommendations of multiple-sample idr
OPTIMALPEAKFILE=$(find "./idr/pairwise/${FACTOR}" -name "mtsc_${FACTOR}*.txt" -type f | xargs wc -l | sort -rn | grep -v ' total$' | head -1 | grep -o './idr.*$')
cp $OPTIMALPEAKFILE idr/optimal/mtsc_$FACTOR\_nexus_idr_peaks.txt
echo "Goodbye!"
