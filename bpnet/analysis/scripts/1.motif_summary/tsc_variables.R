#Melanie Weilert
#January 2021
#Purpose: Define variables and filepaths relevant for analysis

tsc.tasks <- c('tfap2c','tead4','cdx2','gata3','yap1')
#tsc.color.vec<-c('#c22f2f','#007d47','#E69F00','#3600e6','#14a0b5')
tsc.color.vec<-c('#C54B9A','#00B3B4','#B1D03E','#373493','#F47F20')
names(tsc.color.vec)<-tsc.tasks

#Contribution filepaths

tsc.contrib.profile.bws<-lapply(tsc.tasks, function(x){
  paste0('/n/projects/kd2200/publication/bpnet/pred_contrib_bws/', 
         x, '.contrib.profile.bw')
})
names(tsc.contrib.profile.bws)<-tsc.tasks

tsc.contrib.counts.bws<-lapply(tsc.tasks, function(x){
  paste0('/n/projects/kd2200/publication/bpnet/pred_contrib_bws/', 
         x, '.contrib.counts.bw')
})
names(tsc.contrib.counts.bws)<-tsc.tasks

# Actual ChIP-nexus filepaths
tsc.actual.norm.bws<-lapply(tsc.tasks, function(x){
  list(
    pos = paste0('/n/projects/kd2200/publication/bpnet_prep/bw/normalized/mtsc_', tolower(x),'_nexus_filtered_combined_normalized_positive.bw'),
    neg = paste0('/n/projects/kd2200/publication/bpnet_prep/bw/normalized/mtsc_', tolower(x),'_nexus_filtered_combined_normalized_negative.bw')
  )
  
})
names(tsc.actual.norm.bws)<-tsc.tasks

tsc.actual.bws<-lapply(tsc.tasks, function(x){
  list(
    pos = paste0('/n/projects/kd2200/publication/bpnet_prep/bw/mtsc_', tolower(x),'_nexus_filtered_combined_positive.bw'),
    neg = paste0('/n/projects/kd2200/publication/bpnet_prep/bw/mtsc_', tolower(x),'_nexus_filtered_combined_negative.bw')
  )
  
})
names(tsc.actual.bws)<-tsc.tasks

# Predicted ChIP-nexus filepaths
tsc.pred.bws<-lapply(tsc.tasks, function(x){
  list(
    pos = paste0('/n/projects/kd2200/publication/bpnet/pred_contrib_bws/', 
                 x,'.preds.pos.bw'),
    neg = paste0('/n/projects/kd2200/publication/bpnet/pred_contrib_bws/', 
                 x,'.preds.neg.bw')
  )
  
})
names(tsc.pred.bws)<-tsc.tasks
