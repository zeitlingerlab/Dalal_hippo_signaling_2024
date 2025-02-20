"""
Melanie Weilert
Stowers Institute
Purpose: Store ancillary functions to applying perturbations
"""

import os
import sys
import json
import pandas as pd
import numpy as np

#Custom functions
sys.path.insert(0, f'/n/projects/mw2098/shared_code/basepairmodels/scripts')
from predict import predict_basepairmodel


def generate_random_seq(seqlen, weights = [.25, .25, .25, .25]):
    """
    Purpose: Generate a random DNA sequence of a specified length.
    """
    import random
    seq = random.choices(['A','C','G','T'], weights = weights, k=seqlen)
    return(''.join(seq))

def motif_coords(motif, position):
    """
    Kudos: I straight up jacked this from Ziga Avsec's bpnet code bpnet.bpnet.simulate#L19
    Purpose: Given motif (string) and a center position, find the motif boundaries.
    """
    start = position - len(motif) // 2
    end = start + len(motif)
    # print(start, end)
    return start, end

def insert_motif(seq, motif, position):
    """
    Kudos: I straight up jacked this from Ziga Avsec's bpnet code bpnet.bpnet.simulate#L25
    Purpose: Given a sequence, inject a motif centered on a position.
    """
    assert position < len(seq)
    start, end = motif_coords(motif, position)
    new_seq = seq[:start] + motif + seq[end:]
    assert len(new_seq) == len(seq)
    return new_seq

def generate_injected_seq(primary_motif, primary_position,
                          secondary_motif='', secondary_distances=[], seqlen=1000,
                          weights = [.25, .25, .25, .25]):
    """
    Kudos: I straight up jacked this from Ziga Avsec's bpnet code bpnet.bpnet.simulate#L37
    Purpose: Given a sequence, inject a motif centered on a position.
        + If you want to inject 2 motifs to compare distances, use:
            + side_motif to define the second motif sequence
            + side_distances to define the second motif center
    """
    random_seq = generate_random_seq(seqlen = seqlen, weights = weights)
    injected_seq = insert_motif(seq = random_seq, motif = primary_motif, position = primary_position)
    if len(secondary_distances)>1:
        print('Warning! You have entered multiple side distances.',
              'This will inject multiple secondary motifs into the SAME sequence.',
              'If you do not want this, loop this function through single side distances.')
    for d in secondary_distances:
        injected_seq = insert_motif(injected_seq, secondary_motif, d)
    return injected_seq

def predict_injected_seq(model, tasks,
                         primary_motif, primary_position,
                         secondary_motif = '', secondary_position = [],
                         window_around_primary_motif = 50,
                         input_seqlen = 2114,
                         output_seqlen = 1000,
                         trials = 512, correct_for_shoulder_synergy = False,
                         use_cpu = True):
    """
    Purpose:
        + Given a primary motif A, inject that motif and predict.
        + Given a secondary motif B, inject that motif with the primary motif as AB and measure synergy.
    Inputs:
        + model: keras model imported by `load_basepairmodel` already
        + tasks: list of labeled tasks in the order of the output profiles
        + primary/secondary motif: string of sequence to inject
        + primary/secondary position: position within OUTPUT window to inject sequences
        + window_around_primary_motif: measure sum and maximum values in a window around primary motif (bp)
        + seqlen: input sequence length
        + trials: how many repeats to predict random sequences
        + correct_for_shoulder_synergy: see below
        + use_cpu: whether to use CPU to predict results
    Outputs: (pd.df of result values, dict of predicted values directly)
    Application:
        1. Measure motif-motif syngery by assigning primary/secondary motifs as A and B, respectively.
            + `correct_for_shoulder_synergy` can be either False or True.
        2. Measure motif vs null distribution by assigning motif as secondary and null as primary.
            + `correct_for_shoulder_synergy` should be False.
        3. When you are doing analysis on factors with multiple channels per TF,
            you can use the output values to merge them after running this function.
    Math:
        1. Inject both motifs (AB) and primary motif (A) and measure AB/A
        2. Correction for shoulder syngery will take the AB component and correct for affinity of task to B.
            + AB - (B - 0), Where:
                - AB: contains both, primary and secondary motif
                - B : contains only secondary motif
                - 0 : doesn't contain any motif
    """

    #Define partial window
    window_min = primary_position - (window_around_primary_motif // 2)
    window_max = primary_position + (window_around_primary_motif // 2)
    assert window_min>=0 and window_max<=output_seqlen, 'Window range and primary motif injection arent compatible.'

    #Compute where the input positions will be
    input_flank = (input_seqlen - output_seqlen) // 2
    input_primary_position = primary_position + input_flank
    input_secondary_position = secondary_position + input_flank

    #Generate repeated trials of sequence injections
    AB_seqs = [generate_injected_seq(primary_motif = primary_motif, primary_position = input_primary_position,
                                     secondary_motif = secondary_motif,
                                     secondary_distances = [input_secondary_position],
                                     seqlen = input_seqlen) for t in range(trials)]
    A_seqs = [generate_injected_seq(primary_motif = primary_motif, primary_position = input_primary_position,
                                   secondary_motif = '', secondary_distances = [],
                                   seqlen = input_seqlen) for t in range(trials)]
    B_seqs = [generate_injected_seq(primary_motif = secondary_motif, primary_position = input_secondary_position,
                                   secondary_motif = '', secondary_distances = [],
                                   seqlen = input_seqlen) for t in range(trials)]
    null_seqs = [generate_injected_seq(primary_motif = '', primary_position = (input_seqlen // 2),
                                   secondary_motif = '', secondary_distances = [],
                                   seqlen = input_seqlen) for t in range(trials)]

    #Add trials together
    seqs = AB_seqs + A_seqs + B_seqs + null_seqs

    #predict profiles
    preds = predict_basepairmodel(model = model, seqs = seqs, use_cpu = use_cpu)

    #apply counts to profile predictions to result in -> [regions x position x task]
    preds_w_counts = np.array([preds[0][i] * preds[1][i] for i in range(trials*4)])

    #separate different injections and average across trials
    labels = ['AB','A','B','null']
    pred_dict = {labels[r]: np.sum(preds_w_counts[(trials*r):(trials*(r+1))], axis = 0) for r in range(4)}


    #Iterate along different tasks
    results_df = pd.DataFrame()
    for i,t in enumerate(tasks):
        #Define reference and alternative arrays for summary collections.
        ref_profile = pred_dict['A'][:, i]
        if correct_for_shoulder_synergy:
            alt_profile = pred_dict['AB'][:, i] - (pred_dict['B'][:, i] - pred_dict['null'][:, i])

        else:
            alt_profile = pred_dict['AB'][:, i]

        #Compute results across the whole window
        results_dict = {}
        results_dict['all_sum_A'] = ref_profile.sum()
        results_dict['all_sum_AB'] = alt_profile.sum()

        max_idx = np.argmax(ref_profile, axis=0)
        results_dict['all_max_A'] = ref_profile.max()
        results_dict['all_AB_at_max_A'] = float(np.squeeze(alt_profile[max_idx]))

        #Compute results across partial window
        ref_w = ref_profile[window_min:window_max]
        alt_w = alt_profile[window_min:window_max]
        results_dict['window_sum_A'] = ref_w.sum()
        results_dict['window_sum_AB'] = alt_w.sum()

        w_max_idx = np.argmax(ref_w, axis=0)
        results_dict['window_max_A'] = ref_w.max()
        results_dict['window_AB_at_max_A'] = float(np.squeeze(alt_w[w_max_idx]))

        #Record metadata
        results_dict['primary_motif'] = primary_motif
        results_dict['primary_position'] = primary_position
        results_dict['secondary_motif'] = secondary_motif
        results_dict['secondary_position'] = secondary_position
        results_dict['task'] = t

        r_df = pd.DataFrame(results_dict, index=[0])
        results_df = results_df.append(r_df)
    return results_df, pred_dict




def generate_alt_sequences(ref_seq, motifs, motif_unique_col, sequence_window_start_col, sequence_window_end_col,
                           comb_max = None, comb_min = None):
    """
    Purpose: Generate all combinations of mutated sequences based on given sequence and desired mutant coordinates.
    Given:
        ref_seq: np.array with shape [l x 4] of the one-hot encoded reference sequence
        motifs: pd.df of the motifs across the reference sequence
        motif_unique_col: column name designating the unique motif label
        sequence_window_start_col: column name designating where in the input sequence the motifs start
        sequence_window_end_col: column name designating where in the input sequence the motifs end
        comb_max: max number of simultaneous mutations allowed
        comb_min: min number of simultaneous mutations allowed
    Output:
        perturb_seqs_dict_list = trials x -> ref/mutation -> seqlen x 4
    """
    import sys
    sys.path.insert(0, f'/n/projects/mw2098/shared_code/basepairmodels/functions')
    from perturb import generate_random_seq
    from predict import one_hot_encode_sequence
    from itertools import combinations, chain

    motifs['width'] = motifs[sequence_window_end_col] - motifs[sequence_window_start_col]

    #Insert relevant mutant sequences in each combination
    def get_alt_seq_from_combo(combo, ref_seq, motifs, motif_unique_col, sequence_window_start_col, sequence_window_end_col, mut_seqs_dict):
        alt_seq = ref_seq.copy()
        for mut in combo:
            motif = motifs[motifs[motif_unique_col]==mut]
            alt_seq[int(motif[sequence_window_start_col].iloc[0]):int(motif[sequence_window_end_col].iloc[0])] = mut_seqs_dict[mut]
        return alt_seq

    #Mark the max depth of combinations
    c_min_depth = 1
    c_max_depth = motifs.shape[0]
    if(comb_max is not None): c_max_depth = comb_max
    if(comb_min is not None): c_min_depth = comb_min
    
    #Ensure these coordinate positions are integers
    motifs[sequence_window_start_col] = motifs[sequence_window_start_col].astype(int)
    motifs[sequence_window_end_col] = motifs[sequence_window_end_col].astype(int)

    #Determine all combinations of unique motifs
    mut_combos = [list(combinations(motifs[motif_unique_col], d)) for d in range(c_min_depth, c_max_depth+1)]
    mut_combos = list(chain(*mut_combos))

    #Record mutant sequences used in this combination: motif -> motif_len x 4
    mut_seqs_dict = {motifs[motif_unique_col].iloc[idx]:
                          one_hot_encode_sequence(generate_random_seq(seqlen = int(motifs['width'].iloc[idx])))
                     for idx in range(motifs.shape[0])}

    # Generate: combo -> 1000 x 4
    alt_seqs_dict = {'_'.join(list(combo)): get_alt_seq_from_combo(combo = combo,
                                                                   ref_seq = ref_seq,
                                                                   motifs = motifs,
                                                                   motif_unique_col = motif_unique_col,
                                                                   sequence_window_start_col = sequence_window_start_col,
                                                                   sequence_window_end_col = sequence_window_end_col,
                                                                   mut_seqs_dict = mut_seqs_dict)
                     for combo in mut_combos}

    #Add reference sequence
    perturb_seqs_dict = {**{'Reference': ref_seq}, **alt_seqs_dict}

    #Clean up results
    return(perturb_seqs_dict)
