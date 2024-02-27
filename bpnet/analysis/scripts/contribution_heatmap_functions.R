
plot_contributions_heatmap<-function(gr, bw, title = "Contributions", reduce = FALSE, upstream = 50, downstream = 50, reverse = FALSE, return_only_df = FALSE, remove_empty_regions = TRUE,
                                     low.fill.color = "#e8360e", mid.fill.color = "white", high.fill.color = "#1a87bd", low.value = -1, high.value = 1, legend.lab = "Contrib. score", x.lab = "Distance to motif (bp)", for_legend_only = FALSE){
  
  #Get signals across regions
  if(reverse){strand(gr) <- ifelse(strand(gr) == '+', '-', '+')}
  if(reduce){gr <- GenomicRanges::reduce(gr)}
  
  message("Importing data for ", length(gr), " regions...")
  mat <- standard_metapeak_matrix(gr, bw, upstream = upstream, downstream = downstream)
  
  if (remove_empty_regions) {
    message("Removing regions with only zeros (regions with no contributions)...", appendLF = FALSE)
    mat <- mat[apply(mat[,-1], 1, function(x) !all(x==0)),]
    message(" Found ", nrow(mat), " regions.")
  }
  
  df <- as_tibble(mat)
  colnames(df) <- seq(upstream - ncol(mat) + 1, downstream)
  df$rownames <- rownames(mat)
  df$order <- seq(nrow(df))
  df <- tidyr::gather(data = df, position, value, colnames(df)[1]:colnames(df)[ncol(df)-2])
  df$position <- as.numeric(df$position)
  message("Squishing data to the given low and high values.")
  df$value_pos <- scales::oob_squish_any(x = df$value, range = c(0, high.value))
  df$value_pos[df$value_pos == 0] <- NA
  df$value_neg <- scales::oob_squish_any(x = df$value, range = c(low.value, 0))
  df$value_neg[df$value_neg == 0] <- NA
  
  if (return_only_df) {
    return(df)
  } else {
    #Plot matrix
    message("Generating the plot...")
    if (!for_legend_only) {
      p <- ggplot() +
        geom_tile(data=df, aes(x=position, y = order, fill = value_pos))+
        ggtitle(title) +
        {if (is.null(mid.fill.color)) scale_fill_gradient(high = high.fill.color, low = low.fill.color, guide = guide_legend(reverse = FALSE), na.value = "transparent", limits = c(0, high.value), space = "Lab") } +
        {if (!is.null(mid.fill.color)) scale_fill_gradient2(high = high.fill.color, low = low.fill.color, mid = mid.fill.color, guide = guide_legend(reverse = FALSE), na.value = "transparent", midpoint = 0, limits = c(0, high.value), space = "Lab") } +
        scale_x_continuous(expand=c(0,0)) + 
        scale_y_continuous(expand=c(0,0)) +
        theme(axis.line = element_line(colour = "black"),
              axis.text.y = element_text(),
              axis.ticks.y = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(),
              panel.background = element_blank(),
              legend.position="right",
              legend.title = element_text())+
        labs(x = x.lab, y = paste(nrow(mat), "regions")) + 
        labs(fill = legend.lab)
      
      n <- ggplot() +
        geom_tile(data=df, aes(x=position, y = order, fill = value_neg))+
        ggtitle(title) +
        {if (is.null(mid.fill.color)) scale_fill_gradient(high = high.fill.color, low = low.fill.color, guide = guide_legend(reverse = FALSE), na.value = "transparent", limits = c(low.value, 0), space = "Lab") } +
        {if (!is.null(mid.fill.color)) scale_fill_gradient2(high = high.fill.color, low = low.fill.color, mid = mid.fill.color, guide = guide_legend(reverse = FALSE), na.value = "transparent", midpoint = 0, limits = c(low.value, 0), space = "Lab") } +
        scale_x_continuous(expand=c(0,0)) + 
        scale_y_continuous(expand=c(0,0)) +
        theme(axis.line = element_line(colour = "black"),
              axis.text.y = element_text(),
              axis.ticks.y = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(),
              panel.background = element_blank(),
              legend.position="right",
              legend.title = element_text())+
        labs(x = x.lab, y = paste(nrow(mat), "regions")) + 
        labs(fill = legend.lab)
      
      data_p <- ggplot_build(p)
      data_n <- ggplot_build(n)
      
      data_p_filter <- data_p$data[[1]][data_p$data[[1]]$fill != "transparent",]
      data_n_filter <- data_n$data[[1]][data_n$data[[1]]$fill != "transparent",]
      
      new_data <- rbind(data_p_filter, data_n_filter)
      new_data <- new_data[with(new_data, order(x, y)),]
      
      data_p$data[[1]] <- new_data
      
      new_plot <- ggplotify::as.ggplot(ggplot_gtable(data_p))
    } else {
      new_plot <- ggplot() +
        geom_tile(data=df, aes(x=position, y = order, fill = value_neg))+
        ggtitle(title) +
        {if (is.null(mid.fill.color)) scale_fill_gradient(high = high.fill.color, low = low.fill.color, guide = guide_legend(reverse = FALSE), na.value = "transparent", limits = c(low.value, high.value), space = "Lab") } +
        {if (!is.null(mid.fill.color)) scale_fill_gradient2(high = high.fill.color, low = low.fill.color, mid = mid.fill.color, guide = guide_legend(reverse = FALSE), na.value = "transparent", midpoint = 0, limits = c(low.value, high.value), space = "Lab") } +
        scale_x_continuous(expand=c(0,0)) + 
        scale_y_continuous(expand=c(0,0)) +
        guides(fill = guide_colourbar()) +
        theme(axis.line = element_line(colour = "black"),
              axis.text.y = element_text(),
              axis.ticks.y = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(),
              panel.background = element_blank(),
              legend.position="right",
              legend.title = element_text())+
        labs(x = x.lab, y = paste(nrow(mat), "regions")) + 
        labs(fill = legend.lab)
    }
    
    
    return(new_plot)
  }
}



plot_contributions_mutual_table<-function(gr_list, bw_list, title = "Contributions", upstream = 50, downstream = 50, low.fill.color = "red", mid.fill.color = "white", high.fill.color = "blue3", legend.lab = "Contrib. score"){
  #Get signals across regions
  message("Importing data...")
  values <- mclapply(gr_list, function(gr){
    if (is.null(upstream) & is.null(upstream)) {
      upstream = unique(width(gr))
      downstream = 1}
    value <- lapply(bw_list, function(bw) {
      median((standard_metapeak_matrix(gr, bw, upstream = upstream, downstream = downstream)))
    })
  }, mc.cores = length(gr_list))
  
  mat <- data.frame(do.call(cbind, values))
  message("Generating the plot...")
  p <- mat %>%
    as.data.frame() %>%
    rownames_to_column("Task") %>%
    pivot_longer(-c(Task), names_to = "Motif", values_to = "counts") %>% mutate_at('counts', as.numeric) %>%
    ggplot(aes(x=Motif, y=Task, fill=counts)) + 
    geom_tile() +
    ggtitle(title) +
    scale_fill_gradient2(high = high.fill.color, low = low.fill.color, mid = mid.fill.color, guide = guide_legend(reverse = FALSE), na.value = "transparent", midpoint = 0, limits = c(ifelse(min(as.numeric(unlist(mat))) > 0, 0, min(as.numeric(unlist(mat)))), max(as.numeric(unlist(mat)))), space = "Lab") +
    guides(fill = guide_colourbar()) +
    theme(axis.line = element_line(colour = "black"),
          axis.text.y = element_text(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          legend.position="right",
          legend.title = element_text())+
    labs(fill = "Median contrib.")
  
  return(p)
}