#!/bin/bash
#$ -N s_3422
#$ -V -j y
#$ -M bioinfo@stowers.org -m ae -b y
#$ -l mem_free=20G,h_vmem=40G

cd /n/core/Bioinformatics/secondary/Zeitlinger/kd2200/MOLNG-3422.mm10.Ens_102/s_mtsc_ezr_a1_td_single_tdbl

STAR --version > star_version.log

# STAR aligner
STAR --readFilesIn /n/analysis/Zeitlinger/kd2200/MOLNG-3422/AAANK3FM5/n_1_1_ATCACG.fastq.gz \
		--genomeDir /n/analysis/indexes/mm10/annotation/Ens_102/STAR_101bp \
		--runThreadN 4 \
		--outSAMtype BAM SortedByCoordinate \
		--outFileNamePrefix s_mtsc_ezr_a1_td_single_tdbl. \
		--outSAMprimaryFlag OneBestScore \
        --outFilterMultimapNmax 20 \
        --outFilterMismatchNoverLmax 0.1 \
        --outFilterType BySJout \
        --alignSJoverhangMin 8 \
        --alignSJDBoverhangMin 1 \
        --outFilterMismatchNmax 999 \
        --alignIntronMin 20 \
        --alignIntronMax 1000000 \
        --alignMatesGapMax 1000000 \
        --readFilesCommand zcat \
        --limitBAMsortRAM 10000000000 \
        --outSAMattributes NH HI MD AS nM \
        --quantMode TranscriptomeSAM GeneCounts

rsem-calculate-expression --version > rsem_version.log

#RSEM
rsem-calculate-expression --no-bam-output --estimate-rspd  --strandedness reverse \
--bam s_mtsc_ezr_a1_td_single_tdbl.Aligned.toTranscriptome.out.bam /n/analysis/indexes/mm10/annotation/Ens_102/RSEM/mm10.Ens_102.RSEM s_mtsc_ezr_a1_td_single_tdbl.RSEM

# ERCC

bowtie2 -x /n/data1/genomes/indexes/ERCC/bowtie2/ERCC92 -p 4 \
		-U /n/analysis/Zeitlinger/kd2200/MOLNG-3422/AAANK3FM5/n_1_1_ATCACG.fastq.gz \
		-S s_mtsc_ezr_a1_td_single_tdbl.sam 2>&1 | tee -a s_mtsc_ezr_a1_td_single_tdbl_ercc.log

samtools view -bS s_mtsc_ezr_a1_td_single_tdbl.sam > s_mtsc_ezr_a1_td_single_tdbl.bam

samtools sort s_mtsc_ezr_a1_td_single_tdbl.bam -o s_mtsc_ezr_a1_td_single_tdbl.sorted.ercc.bam

samtools index s_mtsc_ezr_a1_td_single_tdbl.sorted.ercc.bam

samtools idxstats s_mtsc_ezr_a1_td_single_tdbl.sorted.ercc.bam > s_mtsc_ezr_a1_td_single_tdbl.ercc.stats

rm s_mtsc_ezr_a1_td_single_tdbl.bam
rm s_mtsc_ezr_a1_td_single_tdbl.sam
rm s_mtsc_ezr_a1_td_single_tdbl.sorted.ercc.bam


Rscript /n/ngs/tools/secundo/scripts/bamTobw_rnaseq.r s_mtsc_ezr_a1_td_single_tdbl.Aligned.sortedByCoord.out.bam reverse



# picard
java -Xmx10g -jar /n/apps/CentOS7/bin/picard.jar CollectRnaSeqMetrics input=s_mtsc_ezr_a1_td_single_tdbl.Aligned.sortedByCoord.out.bam \
  output=s_mtsc_ezr_a1_td_single_tdbl.rnaseq.stats strand_specificity=FIRST_READ_TRANSCRIPTION_STRAND \
  ref_flat=/n/analysis/indexes/mm10/annotation/Ens_102/tables/mm10.Ens_102.refFlat.txt ribosomal_intervals=/n/analysis/indexes/mm10/annotation/Ens_102/extras/mm10.Ens_102.riboList.default.txt 2>&1 | tee -a s_mtsc_ezr_a1_td_single_tdbl_picard.log
