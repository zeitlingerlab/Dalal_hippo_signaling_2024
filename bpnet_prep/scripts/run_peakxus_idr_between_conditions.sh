#!bin/bash

# Melanie Weilert
# October 2019
# Purpose: Given a string representing a factor, run between-condition pairwise (non-repeating) combinations of IDR.

# Notes: The PeakXus must be converted to a BED from a .igv
# Notes: This assumes that there IS an "oracle" peak set to compare to.
# Notes: Make sure your within-condition samples are compariable.
# Notes: Running this script repeatedly can corrupt the .log into making weird statements. Clear idr file environment before rerunning.

# Based on Snakemake file setup.

PEAKXUS_FILE_A=$1 # input sample name
IDR_OUTPUT_DIR=$2 # output idr directory
CONDITION=$3 #prefix string indicating desired condition to grep for
PEAKXUS_COMBINED_FILE=$4 #PEAKXUS combined between-replicates file
PEAKXUS_FILEPATH=$5 #Filepath that contains all of the peakxus bed files (../A.bed)

# Define variables (subset PEAKXUS more strictly than "across conditions")
# PEAKXUS_FILEPATH=$(echo ${PEAKXUS_FILE_A} | sed 's!\(.*\)/.*!\1!')
PEAKXUS_FILES_B=$(ls ${PEAKXUS_FILEPATH}/${CONDITION}*/${CONDITION}*all_transition_points.narrowPeak) #get filepath, then select all that contain condition prefix
PREFIX_A=$(echo ${PEAKXUS_FILE_A} | sed 's!.*/!!' | sed 's/_all_transition_points.*//') #echo, basename, get everything before `all_transition_points`

#Make the IDR if it is not existing
mkdir ${IDR_OUTPUT_DIR}

echo "Running IDR for ${PREFIX_A}..."
#Run pairwise replicates of idr
for PEAKXUS_FILE_B in $PEAKXUS_FILES_B; do
  echo ${PEAKXUS_FILE_A}
  echo ${PEAKXUS_FILE_B}
  PREFIX_B=$(echo ${PEAKXUS_FILE_B} | sed 's!.*/!!' | sed 's/_all_transition_points.*//') #echo, basename, get everything before `all_transition_points`

  #Check for combination
  IDR_OUTPUT_A_B=${PREFIX_A}_vs_${PREFIX_B}_idr.txt
  IDR_OUTPUT_B_A=${PREFIX_B}_vs_${PREFIX_A}_idr.txt

  #Write to log

  echo -e "\n{${IDR_OUTPUT_A_B}}" >> ${IDR_OUTPUT_DIR}/${PREFIX_A}.log

  if [ ${PEAKXUS_FILE_A} != ${PEAKXUS_FILE_B} ]; then #don't run file with itself
    if [[ -f ${IDR_OUTPUT_DIR}/${IDR_OUTPUT_A_B} ]] || [[ -f ${IDR_OUTPUT_DIR}/${IDR_OUTPUT_B_A} ]]; then #if the file exists, do nothing
      echo -e "Skipped because file exists." >> ${IDR_OUTPUT_DIR}/${PREFIX_A}.log
      echo "File exists, skipping IDR."
    else
      echo -e "IDR ran on this combo." >> ${IDR_OUTPUT_DIR}/${PREFIX_A}.log
      idr --samples ${PEAKXUS_FILE_A} ${PEAKXUS_FILE_B} --idr-threshold 0.05 \
      --output-file ${IDR_OUTPUT_DIR}/${IDR_OUTPUT_A_B} \
      --peak-list ${PEAKXUS_COMBINED_FILE}
    fi
  else
    echo -e "Skipped because file is same." >> ${IDR_OUTPUT_DIR}/${PREFIX_A}.log
  fi

done

echo "Goodbye!"
