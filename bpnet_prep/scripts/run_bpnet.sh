#!/bin/bash

#Inputs
DATASPECS=$1 #comma separated list of BPNet dataspecs to go through
GINS=$2 #comma separated list of BPNet gin configurations to go through
MODEL_OUTDIR=$3 #base filepath to write all models to (model)
CONTRIB_OUTDIR=$4 #base filepath to write all contribution scores to (contrib)

for DATASPEC in $(echo $DATASPECS | sed "s/,/ /g")
do
    for GIN in $(echo $GINS | sed "s/,/ /g")
    do
      echo "Running BPNet for data: ${DATASPEC}, configuration: ${GIN}..."
      # OUTDIR="model"
      DATASPEC_ID=$(basename $DATASPEC '_dataspec.yaml') #[dataspec]_[config]
      CONFIG_ID=$(basename $GIN '.gin')
      ID="${DATASPEC_ID}_${CONFIG_ID}"

      #Train model
      #bpnet train --run-id ${ID} --config ${GIN} --memfrac-gpu .9 --vmtouch --in-memory ${DATASPEC} ${MODEL_OUTDIR}
     #Generate contributions
      bpnet contrib --memfrac-gpu .9 ${MODEL_OUTDIR}/${ID} ${CONTRIB_OUTDIR}/${ID}/contrib.h5
      #Generate null contributions
      bpnet contrib --shuffle-seq --memfrac-gpu .9 ${MODEL_OUTDIR}/${ID} ${CONTRIB_OUTDIR}/${ID}/contrib_null.h5

    done
done

echo "Done! =)"
