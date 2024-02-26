---
Date: Feb 2024
Purpose: Give setup directions for Hippo signaling paper
---
Thanks, and Kudos to Melanie Weilert for many scripts!

# Introduction

Here, we give setup instructions for someone interested in reproducing analysis from the Hippo signaling paper. The steps we take are:

1. Environment setup (Anaconda3)
2. Process sequencing data (Snakemake)
3. Analysis (Python and R)

# Conda (Anaconda3) environment setup

BPNet use TensorFlow1 It is HIGHLY recommended that you use `conda=4.7.12` while reproducing this code, since BPNet operates under older sotfware and updated versions of conda are not compatible with the setup instructions below:

## Setup for BPNet environment

The BPNet conda environment can be installed using the instructions found here: [https://github.com/kundajelab/bpnet]. There are 2 environments: 1 with and 1 without a GPU capability. If you choose to install the GPU-compatible BPNet environment on an Nvidia GPU (we trained on a NVIDIA® TITAN RTX GPU), then you will need the appropriate drivers:

+ CUDA v9.0
+ cuDNN v7.0.5

# Process sequencing data

All data is located in `data_preperation/*` and the pipeline instructions are designated from two `Snakefiles` using `Snakemake`. The `Snakefile` sources all the input starting information from the `setup/samples.csv` and `setup/samples_rna.csv` file from the `starting_file` column.

Transcription factor binding datasets were processed separately to train BPNet model with Snakefile under `bpnet_prep/` folder, which I plan to merge with other data.

In order to assign the nexus barcodes, we should parse through each site to get sequencing data.

```
parallel -j 10 bash scripts/nexus_identify_fixed_barcodes.sh -i {} -o txt/nexus_barcodes/\`basename {} .fastq.gz \`\.freqs.txt ::: fastq/mm10/nexus/*.fastq.gz
tail -n +1 txt/nexus_barcodes/*.str.txt
```

In order to process the data, navigate to the `data_preparation/` folder, then type `snakemake -j 6` for 6 simultaneous tasks running. *NOTE: currently broken and will fix asap.

## Software versions associated with data processing

+ R==4.2.0
+ Python==3.7.6
+ bowtie1==1.1.2
+ bowtie2==2.3.5.1
+ cutadapt==2.5
+ samtools==1.14
+ Java OpenJDK==1.8.0_191
+ bamCompare==3.1.3
+ macs2==2.2.6
+ idr==2.0.3
+ snakemake==5.4.5

# Analysis

The rendered .ipynb and .Rmd files are under the `analysis/` folder. Files are numbered in the order by which they were run. Raw figures can be found here as well as code and associated scripts to run analysis.

# The prep and analysis in this folder were done in following order:
1) The `bpnet_prep/` folder is where I had processed binding data with snamekmake and was used to train BPNet model
2) The `data_preperation/` folder is where I ran snakemake to generate all bws of different datasets including binding data.
3) The `bpnet/` folder has BPNet model outputs and analyses associated with it.
