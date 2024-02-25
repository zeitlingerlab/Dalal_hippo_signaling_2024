"""
Melanie Weilert
Stowers Institute
Purpose: Store ancillary functions to modifying/generating plots from python coordinates
"""

import sys
import os
import json
import pandas as pd
import numpy as np
import plotnine
from plotnine import *

def plot_metapeaks(coords_df, input_json,
                   upstream = 250, downstream  = 250,
                   savefig = None, return_df = False,
                   figure_size = (10,6), flip_negative_coverage_strand = True):
    """
    Purpose: Given a pd.df of coordinates (chrom, start, end, strand) and a
        read-in json nested-dictionary (matches basepairmodels.train format),
        generate a metapeak across regions.
    Input:
        + coords_df: pd.df with chrom, start, end, strand columns
        + input_json: directory with given bigwigs to plot and associated grouping information (task -> {bigwig_path, task_id, task_name, strand})
        + upstream: window boundary upstream of the 5' -> 3' oriented coordinate (from the coordinate CENTER)
        + downstream: window boundary downstream of the 5' -> 3' oriented coordinate (from the coordinate CENTER)
        + savefig: if not None , plot the figure to the desired output directory
        + return_df: if not False, return only pd.df of tidied signal values
        + figure_size: dimensions to write image to
        + flip_negative_coverage_strand: if input_json contains multiple strands across a task_id, if True, multiply negative stranded signals by -1
    Output:
        + either a plotnine object or a pd.df of data results
    """
    sys.path.insert(0, f'/n/projects/mw2098/shared_code/basepairmodels/functions')
    from coverage import extract_oriented_coverage
    import plotnine
    plotnine.options.figure_size = figure_size

    #For each task (regardless of strand), generate array across seqlet
    mps_df = pd.DataFrame()
    for t in input_json.keys():
        regions_arr = extract_oriented_coverage(coords_df = coords_df, bigwig_path = input_json[t]['signal'],
                                                upstream = upstream, downstream = downstream)
        mp_vec = np.mean(regions_arr, axis = 0)
        mp_df = pd.DataFrame([mp_vec]).transpose()
        mp_df.columns = ['signal']
        mp_df['position'] = list(range(upstream*-1, downstream))
        mp_df['task_id'] = input_json[t]['task_id']
        mp_df['strand'] = input_json[t]['strand']
        mp_df['task_name'] = input_json[t]['task_name']
        mp_df['task'] = t
        mps_df = mps_df.append(mp_df)

    if flip_negative_coverage_strand:
        mps_df.loc[mps_df['strand']==1, 'signal'] = mps_df.loc[mps_df['strand']==1, 'signal']*-1

    if return_df:
        return(mps_df)
    else:
        mp_plot = (ggplot(data = mps_df, mapping = aes(x = 'position', y = 'signal', group = 'strand'))+
             geom_line()+
             scale_x_continuous(name = 'Position (bp)')+
             scale_y_continuous(name = 'Average signal')+
             facet_grid('. ~ task_name') +
             theme_classic()
                  )
        if savefig:
            outdir = os.path.dirname(savefig)
            os.makedirs(outdir, exist_ok=True)
            mp_plot.save(savefig, verbose = False)
        return(mp_plot)


def plot_heatmaps(coords_df, input_json, order_by_sum = True,
                  removal_threshold = .5, normalize_threshold = .99,
                   upstream = 250, downstream  = 250,
                   savefig = None, return_df = False,
                   figure_size = (10,6)):
    """
    Purpose: Given a pd.df of coordinates (chrom, start, end, strand) and a
        read-in json nested-dictionary (matches basepairmodels.train format),
        generate a normalized heatmap across regions.
    Input:
        + coords_df: pd.df with chrom, start, end, strand columns
        + input_json: directory with given bigwigs to plot and associated grouping information (task -> {bigwig_path, task_id, task_name, strand})
        + order_by_sum: if True [default], order regions by most unnormalized signal across regions
        + removal_threshold: Region quantile to remove all values below
        + normalize_threshold: Region quantile to cap all values above (i.e. applies ceiling to outliers)
        + upstream: window boundary upstream of the 5' -> 3' oriented coordinate (from the coordinate CENTER)
        + downstream: window boundary downstream of the 5' -> 3' oriented coordinate (from the coordinate CENTER)
        + savefig: if not None , plot the figure to the desired output directory
        + return_df: if not False, return only pd.df of tidied signal values
        + figure_size: dimensions to write image to
    Output:
        + either a plotnine object or a pd.df of data results
    """
    sys.path.insert(0, f'/n/projects/mw2098/shared_code/basepairmodels/functions')
    from coverage import extract_oriented_coverage, normalize_coverage_matrix
    import plotnine
    plotnine.options.figure_size = figure_size

    #For each task (regardless of strand), generate array across seqlet
    hms_df = pd.DataFrame()
    for t in input_json.keys():
        #Extract regions in a matrix
        regions_arr = extract_oriented_coverage(coords_df = coords_df, bigwig_path = input_json[t]['signal'],
                                                upstream = upstream, downstream = downstream)
        if order_by_sum:
            regions_arr = regions_arr[np.argsort(np.sum(regions_arr, axis = 1))[::-1], :]

        #Normalized matrix
        regions_norm_arr = normalize_coverage_matrix(regions_arr,
                                                     removal_threshold = removal_threshold,
                                                     normalize_threshold = normalize_threshold)

        #Format matrix into a heatmap pd.df
        hm_df = pd.DataFrame(regions_norm_arr)
        hm_df['row'] = list(range(hm_df.shape[0]))
        hm_melt_df = hm_df.melt(id_vars = ['row'], var_name = 'position', value_name = 'signal')
        hm_melt_df = hm_melt_df.dropna()
        hm_melt_df['position'] = hm_melt_df['position'].astype(int) - upstream

        hm_melt_df['task_id'] = input_json[t]['task_id']
        hm_melt_df['strand'] = input_json[t]['strand']
        hm_melt_df['task_name'] = input_json[t]['task_name']
        hm_melt_df['task'] = t
        hms_df = hms_df.append(hm_melt_df)

    if return_df:
        return(hms_df)
    else:
        #If double stranded, then develop a different color and plotting scale
        if len(set(list(hms_df.strand)))==2:
            colors = ['#08306b', '#08519c' , '#2171b5', '#4292c6', '#6baed6',
                      '#9ecae1', '#c6dbef','#deebf7', '#f7fbff', "white",
                      '#fff5f0','#fee0d2','#fcbba1','#fc9272','#fb6a4a',
                      '#ef3b2c','#cb181d','#a50f15','#67000d']
            rev_strand = sorted(list(set(list(hms_df.strand))))[1]
            # hms_df.signal[hms_df.loc[hms_df.strand==rev_strand, 'signal']] = hms_df.signal[hms_df.loc[hms_df.strand==rev_strand, 'signal']] *-1
            hms_df.loc[hms_df.strand==rev_strand, 'signal'] = hms_df.loc[hms_df.strand==rev_strand, 'signal'] *-1

            #TODO: Test this plot on double stranded data.
            hm_plot = (ggplot(data = hms_df, mapping = aes(x = 'position', y = 'row', fill = 'signal'))+
                geom_tile(alpha = .5)+ #perhaps makes alpha column based on strand?
                scale_fill_gradientn(colors = colors, name = 'Norm. signal')+
                scale_x_continuous(name = 'Position (bp)')+
                scale_y_reverse(name = 'Regions')+
                facet_grid('. ~ task_name') +
                theme_classic())
        else:
            hm_plot = (ggplot(data = hms_df, mapping = aes(x = 'position', y = 'row', fill = 'signal'))+
                geom_tile()+
                scale_fill_gradient(high = '#b2182b', low = 'white', name = 'Norm. signal', limits = [0,1])+
                scale_x_continuous(name = 'Position (bp)')+
                scale_y_reverse(name = 'Regions')+
                facet_grid('. ~ task_name') +
                theme_classic())

        if savefig:
            outdir = os.path.dirname(savefig)
            os.makedirs(outdir, exist_ok=True)
            hm_plot.save(savefig, verbose = False)
        return(hm_plot)


def plot_sequences(coords_df, fasta_path,
                   upstream = 250, downstream  = 250,
                   savefig = None, return_df = False,
                   figure_size = (10,6)):
    """
    Purpose: Given a pd.df of coordinates (chrom, start, end, strand), generate a plot of underlying seuqences.
    Note: These sequences will be oriented in a 5' -> 3' direction relative to the coordinate strand.
    Input:
        + coords_df: pd.df with chrom, start, end, strand columns
        + input_json: directory with given bigwigs to plot and associated grouping information (task -> {bigwig_path, task_id, task_name, strand})
        + fasta_path: path to .fasta file matching genome of motifs
        + upstream: window boundary upstream of the 5' -> 3' oriented coordinate (from the coordinate CENTER)
        + downstream: window boundary downstream of the 5' -> 3' oriented coordinate (from the coordinate CENTER)
        + savefig: if not None , plot the figure to the desired output directory
        + return_df: if not False, return only pd.df of tidied signal values
        + figure_size: dimensions to write image to
        + flip_negative_coverage_strand: if input_json contains multiple strands across a task_id, if True, multiply negative stranded signals by -1
    Output:
        + either a plotnine object or a pd.df of data results
    """
    from Bio.Seq import Seq
    sys.path.insert(0, f'/n/projects/mw2098/shared_code/basepairmodels/functions')
    from motifs import extract_seqs_from_df, resize_coordinates
    import plotnine
    plotnine.options.figure_size = figure_size

    # Resize coordinates to desired window
    coords_resized_df = resize_coordinates(coords_df, width = 1, fix = 'center')
    coords_resized_df = resize_coordinates(coords_resized_df, width = upstream, fix = 'end')
    coords_resized_df = resize_coordinates(coords_resized_df, width = upstream + downstream, fix = 'start').reset_index()

    # Extract sequences on forward orientation of genome
    seqs_unoriented_list = extract_seqs_from_df(coords_resized_df, fasta_path, chrom_column = 'chrom', start_column = 'start', end_column = 'end')

    # Orient sequences
    seqs_oriented_list = [str(Seq(s).reverse_complement()) if coords_resized_df.loc[i, 'strand']=='-' else s for i,s in enumerate(seqs_unoriented_list)]

    seqs_df = pd.DataFrame()
    for i,s in enumerate(seqs_oriented_list):
        s_df = pd.DataFrame(list(s), columns = ['nt'])
        s_df['position'] = list(range(-upstream, downstream))
        s_df['region'] = i
        seqs_df = seqs_df.append(s_df)
    seqs_df['nt'] = pd.Categorical(seqs_df['nt'], ordered=False, categories=['A','C','G','T','N'])

    if return_df:
        return(seqs_df)
    else:
        seq_plot = (ggplot(data = seqs_df, mapping = aes(x='position', y='region', fill='nt'))+
            geom_tile()+
            scale_y_reverse(name = 'Region')+ #order region 1:end from top:bottom
            scale_x_continuous(name = 'Position (bp)')+
            scale_fill_manual(values = ["#36982F", "#402CFD", "#FFB530", "#FC3437", "#808080"], name="Nucleotide")+
            theme_classic()+
            theme(legend_position="right", panel_background = element_blank()))
        if savefig:
            outdir = os.path.dirname(savefig)
            os.makedirs(outdir, exist_ok=True)
            seq_plot.save(savefig, verbose = False)
        return(seq_plot)


def plot_contribution(coords_df, shap_bw_paths_expression, fasta_path,
                   savefig = None,
                   subplot_figure_size = (12,2)):
    """
    Purpose: Given a pd.df of coordinates (chrom, start, end, strand), generate a plot of underlying seuqences.
    Note: These sequences will be oriented in a 5' -> 3' direction relative to the coordinate strand.
    Input:
        + coords_df: pd.df with chrom, start, end, strand columns
        + shap_bw_paths: glob expression to input the different shap .bw files to use for contribution
        + savefig: if not None , plot the figure to the desired output directory
        + subplot_figure_size: dimensions to write image to. Each subplot will be this size.
    Output:
        + either a plotnine object or a pd.df of data results
    """
    import glob
    from Bio import motifs
    from Bio.Seq import Seq
    from predict import one_hot_decode_sequence
    import logomaker
    import matplotlib.pyplot as plt

    sys.path.insert(0, f'/n/projects/mw2098/shared_code/basepairmodels/functions')
    from coverage import extract_oriented_coverage
    from motifs import extract_seqs_from_df
    from predict import one_hot_encode_sequences

    width = list(set(list(coords_df.end - coords_df.start)))[0]
    shap_bw_paths = glob.glob(shap_bw_paths_expression)

    # make Figure and Axes objects
    half_paths = int(len(shap_bw_paths)/2)

    if len(shap_bw_paths)>=2:
        fig, axs = plt.subplots(half_paths, 2, figsize=[subplot_figure_size[0],(half_paths)*subplot_figure_size[1]])
    else:
        fig, axs = plt.subplots(1, 1, figsize=subplot_figure_size)

    contrib_max = []
    for i,bw in enumerate(shap_bw_paths):

        # Extract sequences on forward orientation of genome
        seqs_unoriented_list = extract_seqs_from_df(coords_df, fasta_path = fasta_path, chrom_column = 'chrom', start_column = 'start', end_column = 'end')

        # Orient sequences
        seqs_oriented_list = [str(Seq(s).reverse_complement()) if coords_df.loc[i, 'strand']=='-' else s for i,s in enumerate(seqs_unoriented_list)]

        #1-hot encode sequences
        seqs = one_hot_encode_sequences(seqs_oriented_list)

        #Extract coverage of data
        shap_cov = extract_oriented_coverage(coords_df = coords_df, bigwig_path = bw, upstream = width//2, downstream = width // 2)

        #Multiply 1he seqs by contribution
        contrib_mat = np.array([np.transpose(np.transpose(seqs[i, :, :]) * shap_cov[i]) for i in range(shap_cov.shape[0])])
        cwm_mat = np.nanmean(contrib_mat, axis = 0)

        contrib_max.append(np.max(cwm_mat))
        #Format into data frame
        logo_df = pd.DataFrame(cwm_mat, columns = ['A','C','G','T'])

        logo_df.index.name = 'pos'
        # create Logo object
        if half_paths==1:
            current_axis = axs[int(i % 2)]
        else:
            current_axis = axs[int(np.floor(i/2)), int(i % 2)]

        logo = logomaker.Logo(logo_df,
                              ax=current_axis,
                              show_spines=True)

        current_axis.set_title(bw)
        fig.tight_layout()

    contrib_max = max(contrib_max)

    #Reassign y-axis so all contribution is relative
    for i,bw in enumerate(shap_bw_paths):
        if half_paths==1:
            current_axis = axs[int(i % 2)]
        else:
            current_axis = axs[int(np.floor(i/2)), int(i % 2)]

        current_axis.set_ylim(top=contrib_max)

    if savefig:
        outdir = os.path.dirname(savefig)
        os.makedirs(outdir, exist_ok=True)
        fig.savefig(savefig)
    plt.close(fig)
    return(fig)
