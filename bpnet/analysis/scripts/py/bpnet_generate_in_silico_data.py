#! /home/kd2200/anaconda3/envs/bpnet-gpu/bin/python

##################################################################
# Computational setup
##################################################################
#Packages
import os
import sys
import pandas as pd
import numpy as np
from tqdm import tqdm
from pybedtools import BedTool
from itertools import permutations

# Settings
os.chdir('/n/projects/kd2200/analysis/BPNet/cegkttyz_tsc_bpnet/grid_search/final_cttgy_tsc_model/bpnet/analysis/files/')
pd.set_option('display.max_columns', 100)
figure_filepath = 'figures/4_binding_insilico_perturbs'

# Custom commands
sys.path.insert(0, f'scripts/py')
from bpnet_data_format_functions import myround, myfloor, myceiling

#Pre-existing variables
model_dir = f'/n/projects/kd2200/analysis/BPNet/cegkttyz_tsc_bpnet/grid_search/final_cttgy_tsc_model/bpnet/model/dataspec.yaml_default_fold_5/'

##################################################################
# Identify motifs
##################################################################

motif_seqs = {
    'tead4': 'ACATTCCTG',
    'tfap2c': 'CCCTCAGGC',
    'cdx2': 'GCCATAAA',
    'gata3': 'AGATAAG',
    'junfos': 'ATGAGTCAT',
    'ctcf': 'CCACTAGGGGGCG',
    'elf5': 'CCGGAAG',
    'gata3_double': 'AGATAAGATCT',
    'tead4_double': 'ACATTCCTGGCATTCC'
}

##################################################################
# Define motif pairs
##################################################################

motif_perms = list(permutations(motif_seqs.keys(), 2))
motif_perms.extend([('tead4', 'tead4'), ('cdx2', 'cdx2'), ('tfap2c', 'tfap2c'),('gata3', 'gata3')])

##################################################################
# Import BPNet model
##################################################################

from bpnet.BPNet import BPNetSeqModel
bpn = BPNetSeqModel.from_mdir(model_dir)  # wrap SeqModel to BPNetSeqModel to get `sim_pred` method

##################################################################
# Define helper functions
##################################################################

from bpnet.simulate import generate_sim

#Get results of bpnet.BPNet.sim_pred into a tidy pd.df
def sim_pred_to_df(sim_pred_profiles):
    dfp_ref = pd.DataFrame()
    for k in sim_pred_profiles.keys():
        df = pd.DataFrame(sim_pred_profiles[k]).reset_index()
        df['task'] = os.path.basename(k)
        df.columns = ['position','pos','neg', 'task']
        dfp_ref = dfp_ref.append(df)
    dfp_ref = dfp_ref.melt(id_vars = ['position','task'], var_name = 'strand', value_name = 'pred')
    return dfp_ref

#Get results of bpnet.simulate.generate_sim into a tidy pd.df
def generate_sim_to_df(generate_sim_profiles):
    dfp_alt = pd.DataFrame()
    for d in range(len(generate_sim_profiles)):
        dfp = sim_pred_to_df(generate_sim_profiles[d][1]['profile'])
        dfp['distance'] = generate_sim_profiles[d][0]
        dfp_alt = dfp_alt.append(dfp)
    return dfp_alt

#Get perturbations and convert into tidy.df (dfs = summary, dfp = profiles)
def get_sim(motif_comb, motif_seqs, side_distances = np.arange(505, 650), center_coords = [475, 525]):
    dfs, profiles = generate_sim(bpn, central_motif=motif_seqs[motif_comb[0]], side_motif=motif_seqs[motif_comb[1]],
                       side_distances=side_distances, center_coords=center_coords, contribution=[], correct=True)
    dfs.central_motif = motif_comb[0]
    dfs.side_motif = motif_comb[1]

    #Get reference profiles
    profiles_ref = bpn.sim_pred(motif_seqs[motif_comb[0]], repeat = 256)
    dfp_ref = sim_pred_to_df(profiles_ref)
    dfp_ref['distance'] = 'Reference'

    #Get perturbed profiles
    dfp_alt = generate_sim_to_df(profiles)

    #Combine
    dfp = dfp_ref.append(dfp_alt)

    return dfs, dfp

##################################################################
# Generate distance-based summaries of data
##################################################################

dfs_all = pd.DataFrame()
for perm in motif_perms:
    dfs, _ = get_sim(motif_seqs = motif_seqs, motif_comb=perm, side_distances = np.arange(505, 650))
    dfs_all = dfs_all.append(dfs)

##################################################################
# Save as an intermediate file
##################################################################

dfs_all.to_csv(f'tsv/perturbs/binding/insilico/insilico_summaries.tsv.gz', sep = '\t')
