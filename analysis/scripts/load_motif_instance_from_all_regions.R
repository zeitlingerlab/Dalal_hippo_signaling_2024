# Whis function will load motif_instance_gr from BPnet modico and cwm_scanning outputs


library(testit)

load_motif_instance_from_all_regions <- function(model_path, tf, pattern_string, filter_TE = TRUE){
  file_path <- paste0(model_path,"/",tf,"/instances_by_pattern","/m0_",pattern_string,"_motif-instances-all-regions.tsv.gz")
  motif_instances.df <- read.table(file = file_path, sep = "\t", header = TRUE)
  gr <- motif_instances.df %>%
    #dplyr::rename(chr = example_chrom, start =pattern_start_abs, end = pattern_end_abs) %>%
    #dplyr::select(c(chr,start,end,strand,example_idx, seq_match_cat,contrib_weighted_cat,contrib_weighted,seq_match_p)) %>%#match_weighted_cat, seq_match_cat, contrib_weighted_cat, contrib_weighted_p, match_weighted_p,match_max_task 
    makeGRangesFromDataFrame(keep.extra.columns = T, starts.in.df.are.0based = T, start.field = 'pattern_start_abs', end.field = 'pattern_end_abs', seqnames.field = 'example_chrom')
  gr <- unique(gr)
  # Add unique instance instance id
  gr$instance_id <- paste(tf, 1:length(gr), sep = "_")
  
  if(filter_TE ==TRUE){
    # Filter TE
    rm <- readRDS("/n/projects/mw2098/genomes/mm10/repeatmasker.mm10.gr.rds")
    erv <- rm[grep("ERV", rm$repeat_class)]
    gr <- gr[which(!overlapsAny(gr, erv, ignore.strand=T))]
  }
  assert("The widths of the motif are not equal. Please check code.",(width(gr) %>% unique %>% length)==1)
  gr
}


# # This function load motif instances with certain pattern from cwm_scanned_instances_all_regions

# library(testit)

# load_all_motif_instances_gr <- function(model_path, tf, pattern_string, filter_TE = TRUE){
#   file_path <- paste0(model_path,tf,"/motif-instances-all-regions.tsv.gz")
#   motif_instances.df <- read.table(file = file_path, sep = "\t", header = TRUE)
#   gr <- motif_instances.df %>%
#     dplyr::filter(pattern==pattern_string) %>% 
#     dplyr::rename(chr = example_chrom, start =pattern_start_abs, end = pattern_end_abs) %>%
#     dplyr::select(c(chr,start,end,strand,example_idx)) %>%
#     makeGRangesFromDataFrame(keep.extra.columns = T, starts.in.df.are.0based = T)
#   gr <- unique(gr) # remove duplicated regions because one motif might be mapped by multiple tasks
#   # Add unique instance instance id
#   gr$instance_id <- paste(tf, 1:length(gr), sep = "_")

#   if(filter_TE ==TRUE){
#     # Filter TE
#     rm <- readRDS("/n/projects/mw2098/genomes/mm10/repeatmasker.mm10.gr.rds")
#     erv <- rm[grep("ERV", rm$repeat_class)]
#     gr <- gr[which(!overlapsAny(gr, erv, ignore.strand=T))]
#   }
#   assert("The widths of the motif are not equal. Please check code.",(width(gr) %>% unique %>% length)==1)
#   gr
# }




#load_motif_instance_from_all_regions <- function(model_path, tf, pattern_string, filter_TE = TRUE){
#  file_path <- paste0(model_path,"/",tf,"/instances_by_pattern","/m0_",pattern_string,"_cwm_scanned_instances_all_regions.tsv.gz")


