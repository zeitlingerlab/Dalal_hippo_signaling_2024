suppressPackageStartupMessages(library(optparse, warn.conflicts=F, quietly=T))

option_list <- list(
  make_option(c("-i", "--igv"),
              type="character",
              help="igv file from PeakXus"),
  make_option(c("-n", "--narrowPeak"),
              type="character",
              help="narrowPeak output file name"))

opt <- parse_args(OptionParser(option_list=option_list))

suppressPackageStartupMessages({
  library(Rsamtools)
  library(rtracklayer)
  library(GenomicAlignments)
  library(magrittr)
})

message("Reading peaks: ", opt$file)
stopifnot(file.exists(opt$igv))
peaks.df<-read.table(opt$igv, header = T, sep = "\t")

message("Making GRanges files... ", opt$file)
peaks.np<-makeGRangesFromDataFrame(peaks.df, keep.extra.columns = F, starts.in.df.are.0based = T)

message("Saving file...")
peaks.np.df<-peaks.df[,c(1,2,3,4,7)]
peaks.np.df$strand<-"."
peaks.np.df$signal<-peaks.df$signal
peaks.np.df$p.value<-peaks.df$p.value
peaks.np.df$q.value<-peaks.df$FDR
peaks.np.df$peak<-resize(peaks.np, 1, "center") %>% start

write.table(peaks.np.df, opt$narrowPeak, quote = F, sep = "\t", row.names = F, col.names = F)
