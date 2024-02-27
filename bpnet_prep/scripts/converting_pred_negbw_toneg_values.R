#Standard packages
"https://bioconductor.org/install"
library(rtracklayer) ; library(GenomicRanges); library(magrittr) ; library(Biostrings)
library(ggplot2) ; library(reshape2); library(plyranges); library(Rsamtools); library(data.table)
 library(plyr); library(viridis); library(dplyr)

#KNITR Option

setwd("/n/projects/kd2200/analysis/BPNet/cegkttyz_tsc_bpnet/cegkttyz_tsc_bpnet_morevalchr/tsc_six_factor_valchr_nobias_output/model/bw/"

file.list <- list.files(path = "/n/projects/kd2200/analysis/BPNet/cegkttyz_tsc_bpnet/cegkttyz_tsc_bpnet_morevalchr/tsc_six_factor_valchr_nobias_output/model/bw/",
             pattern = "*.tsc_edit.preds.neg.bw",
             full.names = T)

pred_neg_bw <- function(pred_neg){
     query <- import(pred_neg)
     mcols(query)$score <- query$score*(-1)
     export.bw(query, gsub("tsc_edit.preds.neg.bw", "tsc_edit2.preds.neg.bw", pred_neg, fix = T))
 }
 
 mclapply(file.list, pred_neg_bw, mc.cores = 8)