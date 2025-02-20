#!bin/bash
module load meme
fimo --oc memesuite --skip-matched-sequence --parse-genomic-coord --max-strand --max-stored-scores 10000000 /n/projects/kd2200/publication/rebuttal/mtsc_tead4_macs2_peaks_nexus_top1000_101bpwidth/meme.txt memesuite/all_combined_tead4_peaks_seqs.fa > memesuite/fimo_matches_default_settings_at_all_combined_tead4_peaks.tsv
