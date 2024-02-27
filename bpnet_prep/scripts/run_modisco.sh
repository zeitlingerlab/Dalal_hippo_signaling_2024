

# Run BPNet
# /n/projects/kd2200/analysis/BPNet/bpnet_processing_streamlined/osknsps_bpnet$ scp -i /n/projects/kd2200/aws/keys/bpnet_updated.pem -r * ec2-user@ec2ec2-3-81-150-238.compute-1.amazonaws.com:workspace/
# wandb login c5fd075d26ba024f5d7a164308eaaafd85d95d08
# bpnet train --wandb-project kd2200/osknsps_bpnet --num-workers 4 --vmtouch --memfrac-gpu .6 --overwrite --run-id osknsps_bpnet dataspec.yaml model_dir (--note-params note='test and vld chr same[2,9,18] in bpnet9.gin')
# bpnet contrib model_dir/osknsps_bpnet/  --method=deeplift --memfrac 0.8 --contrib-wildcard='*/profile/wn' model_dir/contrib_file/imp_score_deeplift.h5
# bpnet contrib model_dir/osknsps_bpnet/  --method=deeplift --memfrac 0.8 --contrib-wildcard='*/profile/wn' --shuffle-seq --max-regions 5000  model_dir/contrib_null_file/imp_score_deeplift_null.h5

#!/bin/bash

TASKS=$1 #comma separated list of tasks to run through

for TASK in $(echo $TASKS | sed "s/,/ /g")
do
   echo "Running TF-MoDISco on $TASK..."
    bpnet modisco-run  model_dir/contrib_file/imp_score_deeplift.h5 --null-contrib-file model_dir/contrib_null_file/imp_score_deeplift_null.h5 --contrib-wildcard=$TASK\/profile/wn \
    --only-task-regions --overwrite  modisco_dir/$TASK\/

    echo "Running ChIP-nexus analysis on $TASK..."
    bpnet chip-nexus-analysis modisco_dir/$TASK --footprint-width=250

    echo "Running CWM-scan for $TASK..."
    bpnet cwm-scan --num-workers 32 --add-profile-features modisco_dir/$TASK modisco_dir/$TASK\/cwm_scanned_seqlet_instances.tsv.gz
done

 # echo "Running CWM-scan for $TASK..."
 #    bpnet cwm-scan modisco_dir/$TASK modisco_dir/$TASK/motif-instances-task-regions.tsv.gz --add-profile-features #`example_idx` is NOT the same across task runs
 #    bpnet cwm-scan --contrib-file model_dir/contrib_file/imp_score_deeplift.h5 modisco_dir/$TASK modisco_dir/$TASK/motif-instances-all-regions.tsv.gz --add-profile-features #`example_idx` is the same across task runs

