"""
Melanie Weilert
Stowers Institute
Purpose: Given a set of inputs, generate metaplots from TF-MoDISco
"""

############################################################################################################
#Setup
############################################################################################################

import os
import json
import sys
import glob
import pandas as pd
import numpy as np
import time
import argparse
import multiprocessing as mp
from tqdm import tqdm
from wand.image import Image as WImage
from pathos.multiprocessing import ProcessingPool as Pool

#Custom functions
sys.path.insert(0, f'/n/projects/mw2098/shared_code/basepairmodels/functions')
from motifs import import_modisco_seqlets_txt, map_modisco_coordinates_to_genome

sys.path.insert(0, f'/n/projects/mw2098/shared_code/bpnet/scripts')
from plots import plot_heatmaps, plot_metapeaks, plot_sequences, plot_contribution

############################################################################################################
#Functions
############################################################################################################

def generate_contrib(seqlet_path, genome, shap_path_expression):
    s_df = pd.read_csv(seqlet_path, sep='\t', header=None, names=['chrom', 'start', 'end', 'name', 'score', 'strand'])
    s_base = os.path.basename(seqlet_path).replace('_seqlets.bed', '')

    #Plot contribution scores
    contrib_output_plot_path = os.path.dirname(seqlet_path) + '/plots/' + s_base + '_contrib.pdf'
    contrib_plot = plot_contribution(coords_df = s_df, shap_bw_paths_expression = shap_path_expression,
                                     fasta_path = genome,
                                     savefig = contrib_output_plot_path)
    contrib_png = WImage(filename=contrib_output_plot_path, resolution=300) # bigger
    contrib_png.format = 'png'
    contrib_png.save(filename=contrib_output_plot_path.replace('.pdf','.png'))
    return(None)

def generate_metapeaks(seqlet_path, bigwig_path, upstream=250, downstream=250, figwidth = 15, figheight = 3):
    s_df = pd.read_csv(seqlet_path, sep='\t', header=None, names=['chrom', 'start', 'end', 'name', 'score', 'strand'])
    s_base = os.path.basename(seqlet_path).replace('_seqlets.bed', '')

    #Generate metapeaks
    mp_output_plot_path = os.path.dirname(seqlet_path) + '/plots/' + s_base + '_metapeaks.pdf'
    mp_plot = plot_metapeaks(coords_df = s_df, bigwig_path = bigwig_path,
                             upstream = upstream, downstream = downstream,
                             savefig = mp_output_plot_path,
                             figure_size = (figwidth,figheight))
    mp_png = WImage(filename=mp_output_plot_path, resolution=300) # bigger
    mp_png.format = 'png'
    mp_png.save(filename=mp_output_plot_path.replace('.pdf','.png'))
    return(None)

def generate_heatmaps(seqlet_path, bigwig_path, upstream=250, downstream=250, figwidth = 15, figheight = 3):
    s_df = pd.read_csv(seqlet_path, sep='\t', header=None, names=['chrom', 'start', 'end', 'name', 'score', 'strand'])
    s_base = os.path.basename(seqlet_path).replace('_seqlets.bed', '')

    #Generate normalized heatmaps
    hm_output_plot_path = os.path.dirname(seqlet_path) + '/plots/' + s_base + '_heatmaps.pdf'
    hm_plot = plot_heatmaps(coords_df = s_df, bigwig_path = bigwig_path,
                            upstream = upstream, downstream = downstream,
                            order_by_sum = True, savefig = hm_output_plot_path,
                            figure_size = (figwidth,figheight))
    hm_png = WImage(filename=hm_output_plot_path, resolution=300) # bigger
    hm_png.format = 'png'
    hm_png.save(filename=hm_output_plot_path.replace('.pdf','.png'))
    return(None)

def generate_seqplots(seqlet_path, genome, bigwig_path, upstream=250, downstream=250):
    import pyBigWig
    s_df = pd.read_csv(seqlet_path, sep='\t', header=None, names=['chrom', 'start', 'end', 'name', 'score', 'strand'])
    s_base = os.path.basename(seqlet_path).replace('_seqlets.bed', '')

    #Collect sequence signal given .bw file
    b = pyBigWig.open(bigwig_path)
    sig = np.stack([b.stats(row.chrom, row.start, row.end) for i,row in s_df.iterrows()])
    bw_arr = np.sum(np.array(sig), axis = 1)
    # print(bw_arr.shape)
    s_df['signal'] = bw_arr

    #Order sequences
    s_df = s_df.sort_values(['signal'], ascending = False)

    seq_output_plot_path = os.path.dirname(seqlet_path) + '/plots/' + s_base + '_seqs.pdf'
    seq_plot = plot_sequences(coords_df = s_df, fasta_path = genome,
                              savefig = seq_output_plot_path,
                              upstream = upstream, downstream = downstream,
                              figure_size = (6,3))
    seq_png = WImage(filename=seq_output_plot_path, resolution=300) # bigger
    seq_png.format = 'png'
    seq_png.save(filename=seq_output_plot_path.replace('.pdf','.png'))
    return(None)

############################################################################################################
#Main code
############################################################################################################

def main():
    # parse the command line arguments
    parser = argparse.ArgumentParser()

    #Add new transfer learning inputs
    parser.add_argument('--seqlet-path-expression', type=str,
                        help="Glob expression to seqlets (usually is '*/*/pattern*_seqlets.bed')")
    parser.add_argument('--bigwig-path', type=str,
                        help="Input .bw filepath required for coverage plots.")
    parser.add_argument('--shap-path-expression', type=str,
                        help="Glob expression to shap .bw contributions (usually is 'f'shap/*/*_scores.bw'')")
    parser.add_argument('--genome', type=str,
                        help="Path to reference genome .fa")
    parser.add_argument('--upstream', type=int,
                        help="How many bp upstream of each motifto plot metapeaks and heatmaps?", default=250)
    parser.add_argument('--downstream', type=int,
                        help="How many bp downstream of each motif to plot metapeaks and heatmaps?", default=250)
    parser.add_argument('--threads', type=int,
                        help="Number of threads to parallelize through", default=16)
    args = parser.parse_args()

    #Import paths to information
    seqlet_paths = sorted(glob.glob(args.seqlet_path_expression))
    seqlet_path_n = len(seqlet_paths)

    #Identify task counts
    tasks = ['atac'] #just one because ChromBPNet doesn't support multi-task
    tasks_unique = list(set([t.replace('_pos', '').replace('_neg', '') for t in tasks]))
    task_unique_count = len(tasks_unique)
    figwidth = task_unique_count*3 #for heatmaps and metapeaks

    #Generate different plots in //


    with Pool(args.threads) as p:
        start = time.time()
        filler = p.map(generate_seqplots, seqlet_paths, [args.genome]*seqlet_path_n, [args.bigwig_path]*seqlet_path_n, [50]*seqlet_path_n, [50]*seqlet_path_n)
        print("Complete")
        end = time.time()
        print('Generating seqplots took (s)= ' + str(end-start))


    # with Pool(args.threads) as p:
    #     start = time.time()
    #     filler = p.map(generate_contrib, seqlet_paths, [args.genome]*seqlet_path_n, [args.shap_path_expression]*seqlet_path_n)
    #     print("Complete")
    #     end = time.time()
    #     print('Generating multi-task contrib logos took (s)= ' + str(end-start))
    #
    # with Pool(args.threads) as p:
    #     start = time.time()
    #     filler = p.map(generate_metapeaks, seqlet_paths, [args.bigwig_path]*seqlet_path_n, [args.upstream]*seqlet_path_n, [args.downstream]*seqlet_path_n, [figwidth]*seqlet_path_n, [3]*seqlet_path_n)
    #     print("Complete")
    #     end = time.time()
    #     print('Generating metapeaks took (s)= ' + str(end-start))
    #
    # with Pool(args.threads) as p:
    #     start = time.time()
    #     filler = p.map(generate_heatmaps, seqlet_paths, [args.bigwig_path]*seqlet_path_n, [args.upstream]*seqlet_path_n, [args.downstream]*seqlet_path_n, [figwidth]*seqlet_path_n, [3]*seqlet_path_n)
    #     print("Complete")
    #     end = time.time()
    #     print('Generating heatmaps took (s)= ' + str(end-start))

if __name__ == '__main__':
    main()
