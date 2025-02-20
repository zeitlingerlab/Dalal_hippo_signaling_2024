"""
Melanie Weilert
Stowers Institute
Purpose: Store ancillary functions to modifying/generating motifs not covered by `basepairmodels.CLI`
"""

import os
import json
import pandas as pd
import numpy as np

def extract_oriented_coverage(coords_df, bigwig_path, upstream = 250, downstream = 250):
    """
    Purpose: Extract a [region x window] array
         + Regions oriented in 5' -> 3' direction
         + Relative to the CENTER of given coordinates.
    Input:
        + coords_df: pd.df with chrom, start, end, strand columns
        + bigwig_path: path to bigwig file
    Output: np.array with [region x window] dimensions and region order matches coords_df
    """
    import pyBigWig
    bw = pyBigWig.open(bigwig_path)

    #Assert that coordinates are the same length
    coord_width = set(list(coords_df.end - coords_df.start))
    assert len(coord_width)==1, 'Input seqlets are different lengths--check input coordinates.'
    coord_half = list(coord_width)[0]//2


    #Define window edges relative to given coordinates
    upstream_flank = upstream - coord_half
    downstream_flank = downstream - coord_half

    coverage_list = []
    for i,row in coords_df.iterrows():

        #Determine strand information
        if row.strand == '-':
            cov = np.flip(bw.values(row.chrom, row.start - downstream_flank, row.end + upstream_flank))
        else:
            cov = bw.values(row.chrom, row.start - upstream_flank, row.end + downstream_flank)
        coverage_list.append(cov)

    return(np.array(coverage_list))


def normalize_coverage_matrix(mat, removal_threshold = .5, normalize_threshold = .99):
    """
    Purpose: Normalize a coverage matrix to quantile minimums and maximums for plotting.
    """
    #Find threshold for each region to cap high values
    quantile_max = np.quantile(mat, normalize_threshold, axis = 1)

    #Find threshold for each region to zero out low values
    quantile_min = np.quantile(mat, removal_threshold, axis = 1)

    mat_norm = mat
    for row in range(mat.shape[0]):
        # Remove regions that do not meet quantile minimum
        mat_norm[row][np.where(mat_norm[row]<=quantile_min[row])] = None

        #Normalize all other values to quantile maximum
        mat_norm[row] = mat_norm[row] / quantile_max[row]

        #Remove outliers by capping higher values than quantile maximum
        mat_norm[row][np.where(mat_norm[row]>1)] = 1

    return(mat_norm)
