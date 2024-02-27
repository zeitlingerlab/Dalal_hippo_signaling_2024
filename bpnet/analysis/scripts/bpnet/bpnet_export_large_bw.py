
"""
Melanie Weilert
February 2020
Purpose: Export prediction to a bigWig file.

This code is modified from github.com/kundajelab/bpnet to handle exporting
large regions alongside GPUs with 16Gb or less memory.

Intructions:
1. Activate `bpnet` conda environment.
2. `python /n/projects/mw2098/shared_code/bpnet/bpnet_export_large_bw.py [options]`
3. For help: `python /n/projects/mw2098/shared_code/bpnet/bpnet_export_large_bw.py -h`
"""

# Setup
import os
import pyBigWig
import logging
import subprocess
import pandas as pd
import numpy as np
import tensorflow as tf
from optparse import OptionParser
from itertools import compress
from operator import itemgetter
from collections import OrderedDict
from tqdm import tqdm
from pybedtools import BedTool
from bpnet.BPNet import BPNetSeqModel
from bpnet.preproc import resize_interval
from bpnet.seqmodel import SeqModel
from bpnet.BPNet import BPNetSeqModel
from bpnet.utils import add_file_logging, read_pkl, create_tf_session
from bpnet.dataspecs import DataSpec
from bpnet.preproc import resize_interval
from argh.decorators import named, arg

import warnings
warnings.filterwarnings("ignore")

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


parser = OptionParser()
parser.add_option("-m", "--model_dir",
                  help="Path to the trained model directory (specified in `bpnet train <output_dir>`")
parser.add_option("-o", "--output_prefix",
                  help="Creates a directory and prefix for BW files. If True, the preexisting output directory will be overwritten")
parser.add_option("-c", "--chrom_sizes",
                  help="Chromosize sizes .txt file by which to merge large bigwig files.")
parser.add_option("-f", "--fasta_file", default = None,
                  help="Fasta file to pull sequences from. If not specified, file specified in dataspec.yml will be used [default = %default]")
parser.add_option("-r", "--regions", default = None,
                  help="Path to the interval bed file. If not specified, files specified in dataspec.yml will be used [default = %default]")
parser.add_option("-d", "--contrib-method", default = 'deeplift',
                  help="Contribution score method to use. Available: grad, deeplift [default = %default]")
parser.add_option("-w", "--contrib-wildcard", default = '*/profile/wn,*/counts/pre-act',
                  help="Wildcard of the contribution scores to compute. For example, */profile/wn computes\n the profile contribution scores for all the tasks (*) using the wn normalization (see bpnet.heads.py).\n*/counts/pre-act computes the total count contribution scores for all tasks w.r.t. the pre-activation output \nof prediction heads. Multiple wildcards can be by comma-separating them. [default = %default]")
parser.add_option("-s", "--scale-contribution", default = False,
                  help="If True, multiple the contribution scores by the predicted count value [default = %default]")
parser.add_option("-b", "--batch-size", default = 512, type = "int",
                  help="Batch size for computing the predictions and contribution scores [default = %default]")
parser.add_option("-g", "--gpu", default = 0, type = "int",
                  help="which gpu to use.  [default = %default]")
parser.add_option("-x", "--memfrac-gpu", default = .45, type = "float",
                  help="what fraction of the GPU memory to use [default = %default]")

(options, args) = parser.parse_args()

def bpnet_export_large_bw(model_dir,
                          output_prefix,
                          chrom_sizes,
                          fasta_file=None,
                          regions=None,
                          contrib_method='grad',
                          contrib_wildcard='*/profile/wn,*/counts/pre-act',  # specifies which contrib. scores to compute
                          batch_size=256,
                          scale_contribution=False,
                          gpu=0,
                          memfrac_gpu=0.45):
    """Export model predictions and contribution scores to big-wig files
    """
    from pybedtools import BedTool
    from bpnet.modisco.core import Seqlet

    regions_file = regions
    # Make directories
    output_dir = os.path.dirname(output_prefix)
    add_file_logging(output_dir, logger, 'bpnet-export-large-bw')
    os.makedirs(output_dir, exist_ok=True)
    if gpu is not None:
        create_tf_session(gpu, per_process_gpu_memory_fraction=memfrac_gpu)

    logger.info("Load model")
    bp = BPNetSeqModel.from_mdir(model_dir)

    #Get list of intervals from either the model_dir or from the provided BED file
    if regions is not None:
        logger.info(f"Computing predictions and contribution scores for provided regions: {regions}")
        regions = list(BedTool(regions))
        regions_df = BedTool(regions_file).sort().to_dataframe() #parse through pd.df for chromosomes
    else:
        logger.info("--regions not provided. Using regions from dataspec.yml")
        ds = DataSpec.load(os.path.join(model_dir, 'dataspec.yml'))
        regions = ds.get_all_regions()
        regions_df = BedTool(regions_file).sort().to_dataframe() #parse through pd.df for chromosomes

    seqlen = bp.input_seqlen()
    logger.info(f"Resizing regions (fix=center) to model's input width of: {seqlen}")
    regions = [resize_interval(interval, seqlen) for interval in regions]

    logger.info("Sort the bed file")
    regions = list(BedTool(regions).sort())

    for chr in regions_df.chrom.unique():
        chr_idx = list(regions_df.chrom == chr)
        chr_idx_for_list = list(compress(range(len(chr_idx)), chr_idx))
        regions_across_chr = list(itemgetter(*chr_idx_for_list)(regions)) #subset list by chromosome
        print("Predicting " + str(len(regions_across_chr)) + " regions across " + chr + "..." )
        bp.export_bw(regions=regions_across_chr,
                    output_prefix=output_prefix + "_" + chr,
                    contrib_method=contrib_method,
                    fasta_file=fasta_file,
                    flip_negative_strand = True,
                    pred_summaries=contrib_wildcard.replace("*/", "").split(","),
                    batch_size=batch_size,
                    scale_contribution=scale_contribution,
                    chromosomes=None)  # infer chromosomes from the fasta file
    logger.info("Chromosome-separated bws generated. Please refer to merge_bigwigs.r for further merging.")
    return None

bpnet_export_large_bw(model_dir = options.model_dir,
                      output_prefix = options.output_prefix,
                      chrom_sizes = options.chrom_sizes,
                      fasta_file = options.fasta_file,
                      regions = options.regions,
                      contrib_method = options.contrib_method,
                      contrib_wildcard = options.contrib_wildcard,
                      batch_size = options.batch_size,
                      scale_contribution = options.scale_contribution,
                      gpu = options.gpu,
                      memfrac_gpu = options.memfrac_gpu)
