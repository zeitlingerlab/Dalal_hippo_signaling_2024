#!bin/bash
module load meme
fimo --oc memesuite --skip-matched-sequence --parse-genomic-coord --max-strand --max-stored-scores 10000000 ./rebuttal/analysis/memesuite/custom_meme.txt ./rebuttal/analysis/memesuite/all_peaks_motifs_ov_seqs.fa > ./rebuttal/analysis/memesuite/fimo_matches_seqlet_pwm.tsv
