#Melanie Weilert
#August 2019
#Purpose: Script to convert IDR to BED files, with formatting options

suppressPackageStartupMessages(library(optparse, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(GenomicRanges, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(BSgenome, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(dplyr, warn.conflicts=F, quietly=T))
suppressPackageStartupMessages(library(rtracklayer, warn.conflicts=F, quietly=T))
source("/n/projects/mw2098/shared_code/rscripts/granges_common.r")

option_list <- list(
  make_option(c("-i", "--idr"), 
              type="character",
              default=NA,
              help="Path of .txt file to process"),
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

#Read IDR file into GRanges
regions.df<-read.table(opt$idr)
colnames(regions.df)<-c("chr", "start", "end", "peak_name", 
                        "scorefordisplay","strand","signal", "pval", "qval", "summit")
regions.gr<-makeGRangesFromDataFrame(df=regions.df, keep.extra.columns = T, starts.in.df.are.0based = T)

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

#Export regions to bed
regions_standard.gr$score<-regions_standard.gr$signal #record average signal across the region as the score column
regions_standard.gr$strand<-"*"

regions_standard.df<-as.data.frame(regions_standard.gr)
write.table(x = regions_standard.df[,1:3], file = opt$output, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
#(regions_standard.gr, opt$output, format="BED")
