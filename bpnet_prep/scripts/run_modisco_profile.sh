#!/bin/bash

#Inputs
DATASPECS=$1 #comma separated list of BPNet dataspecs to go through
GINS=$2 #comma separated list of BPNet gin configurations to go through
MODEL_OUTDIR=$3 #base filepath to write all models to (model)
CONTRIB_OUTDIR=$4 #base filepath to write all contribution scores to (contrib)
TASKS=$5 #comma separated list of BPNEt tasks to loop through
MODISCO_OUTDIR=$6  #base filepath to write all modisco scores to (modisco)

for DATASPEC in $(echo $DATASPECS | sed "s/,/ /g")
do
    for GIN in $(echo $GINS | sed "s/,/ /g")
    do
        for TASK in $(echo $TASKS | sed "s/,/ /g")
        do

          echo "Running MODISCO for data: ${DATASPEC}, configuration: ${GIN}, task: ${TASK}..."
          # OUTDIR="model"
          DATASPEC_ID=$(basename $DATASPEC '_dataspec.yaml') #[dataspec]_[config]
          CONFIG_ID=$(basename $GIN '.gin')
          ID="${DATASPEC_ID}_${CONFIG_ID}"

          #Get modisco
          bpnet modisco-run --null-contrib-file ${CONTRIB_OUTDIR}/${ID}/contrib_null.h5 \
          --contrib-wildcard=${TASK}\/profile/wn --only-task-regions \
          ${CONTRIB_OUTDIR}/${ID}/contrib.h5 ${MODISCO_OUTDIR}/${ID}/${TASK}

        #   echo "Running ChIP-nexus analysis on $TASK..."
        #   bpnet chip-nexus-analysis ${MODISCO_OUTDIR}/${ID}/${TASK} --footprint-width=800

        #   echo "Running CWM-scan for $TASK..."
        #   bpnet cwm-scan --contrib-file ${CONTRIB_OUTDIR}/${ID}/contrib.h5 ${MODISCO_OUTDIR}/${ID}/${TASK} \
        #   ${MODISCO_OUTDIR}/${ID}/${TASK}/motif-instances-all-regions.tsv.gz #`example_idx` is the same across task runs
        #   bpnet cwm-scan ${MODISCO_OUTDIR}/${ID}/${TASK} \
        #   ${MODISCO_OUTDIR}/${ID}/${TASK}/motif-instances-task-regions.tsv.gz #`example_idx` is NOT the same across task runs

         done
    done
done

echo "Done! =)"
