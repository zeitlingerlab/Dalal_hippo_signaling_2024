

#Khyati 
#September 2019
#Purpose: Script to combine and reduce IDR replicates, with formatting options
suppressPackageStartupMessages(library(optparse, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(GenomicAlignments, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(GenomicRanges, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(utils, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(BSgenome, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(dplyr, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(rtracklayer, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(data.table, warn.conflicts=F, quietly=T))
source("/n/projects/mw2098/shared_code/rscripts/granges_common.r")


# idr_reduce_combine.r - some commands for processing replicate idr

option_list <- list(
  make_option(c("-i", "--idr"), 
              type="character",
              default=NA,
              help="Path of idr.txt files to process"),
  make_option(c("-b", "--blacklist"), 
              type="character",
              default=NA,
              help="Filepath to the blacklisted genome"),
  make_option(c("-g", "--genome"), 
              type="character",
              default=NA,
              help="Genome of sample (correct string can be found by running `names(genomeStyles())` in R.)"),
  make_option(c("-w", "--boundaryWidth"), 
              type="integer",
              default=1200,
              help="Padding allowance when checking peak proximity to chromosome ends"),  
  make_option(c("-o", "--output"),
              type="character",
              default=NA,
              help="Output BED filepath")
)
# ------------------------------------------------------------------------------
#Get options
opt <- parse_args(OptionParser(option_list=option_list))


#Read IDR.txt files for replicates
#regions.df =  list.files(path= opt$idr, pattern= "*idr.txt")  
regions.df <- list.files(path= opt$idr, pattern="*.txt", full.names=TRUE)
region_df <- lapply(regions.df, read.table)
#regions.df <- lapply(regions.df, function(x) {read.table(file = x, header = FALSE)})
#combined_df <- do.call("rbind", lapply(regions.df))
combined_df <- rbindlist(region_df, use.names=FALSE)
colnames(combined_df)<-c("chr", "start", "end", "string", "idr_score", 
                        "strand", "signal", "pval", "qval", "summit", 
                        "local_idr", "global_idr", "rep1_start", "rep1_end", 
                        "rep1_signal", "rep1_summit", "rep2_start", "rep2_end", 
                        "rep2_signal", "rep2_summit")
regions.gr<-makeGRangesFromDataFrame(df=combined_df, keep.extra.columns = T, starts.in.df.are.0based = T)

#Blacklist regions
blacklist.gr<-import(opt$blacklist)
ov<-findOverlaps(regions.gr, blacklist.gr, ignore.strand=T)
if(length(ov)>0){regions.gr<-regions.gr[-ov@from]}

#Assign genome
genome<-getBSgenome(opt$genome)
seqlevels(regions.gr)<-seqlevels(genome)
seqinfo(regions.gr)<-seqinfo(genome)

#Check chromosome boundaries
regions_bounded.gr<-check_chromosome_boundaries(regions.gr, resize_boundary = opt$boundaryWidth, genome=genome)

#Keep standard chromosomes
regions_standard.gr<-keepStandardChromosomes(regions_bounded.gr, pruning.mode = "coarse")

#Format regions to bed
regions_reduced.gr <- GenomicRanges::reduce(regions_standard.gr, ignore.strand= TRUE) %>% resize(., width=1, fix="center")
regions_reduced.gr$strand<-"*"

# export as bed
regions_reduced.gr<-as.data.frame(regions_reduced.gr)
write.table(x = regions_reduced.gr[,1:3], file = opt$output, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
#export.bed(regions_reduced.gr, opt$output)



