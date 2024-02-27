# Generating unstranded polII ChIP-nexus bigwig file

gr  <- readRDS("/n/projects/kd2200/analysis/BPNet/cegkttyz_tsc_bpnet/grid_search/cttgy_tsc_model_macs2_biastrack_valchrom_v1/rdata/nexus/mtsc_polii_nexus_filtered_combined.granges.rds")
gr <- resize(gr, 1, "start")
gr.cov <- coverage(gr)
export(gr.cov,"/n/projects/kd2200/analysis/BPNet/cegkttyz_tsc_bpnet/grid_search/cttgy_tsc_model_macs2_biastrack_valchrom_v1/bw/mtsc_polii_nexus_combined_unstranded.bw")