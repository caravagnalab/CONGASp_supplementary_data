#VECCHIO
# With respect to the original function, this one plots only the histogram of clusters
# identified by CONGAS+
plot_fit_density_mod = function (x, highlights = TRUE, position = 'stack', palette_name = 'Set1') 
{
  CNAs = get_fit(x, "CNA")
  if (highlights) {
    nclusters = CNAs$cluster %>% unique %>% length()
    CNAs = CNAs %>% group_by(segment_id, value) %>% mutate(grp_size = n()) %>% 
      filter(grp_size != nclusters) %>% pull(segment_id) %>% 
      unique
    cli::cli_alert("Plotting segments where different CNAs are present: {.field {CNAs}}.")
  }
  else {
    CNAs = CNAs %>% pull(segment_id) %>% unique
    cli::cli_alert("Showing all segments (this plot can be large).")
  }
  if (length(CNAs) == 0) {
    return(ggplot() + geom_blank())
  }
  plots <- Rcongas:::plot_data_histogram(x, CNAs)
  data_hist <- plots$data
  clustering_assignments = get_fit(x, what = "cluster_assignments") %>% 
    dplyr::select(-modality)
  data_hist$modality <- sapply(data_hist$modality %>% strsplit(" "), 
                               function(y) y[1])
  data_hist <- dplyr::left_join(data_hist, clustering_assignments)
  mixing_props <- Rcongas:::get_mixing_proportions(x)
  densities <- Rcongas:::assemble_likelihood_tibble(x, CNAs)
  colnames(densities)[c(3, 5)] <- c("cluster", "segment_id")
  densities <- dplyr::left_join(densities, mixing_props) %>% 
    mutate(value = value * (mixing + 0.001))
  CNA <- Rcongas:::get_CNA(x)
  colnames(CNA)[3] <- "CNA"
  data_hist <- dplyr::left_join(data_hist, CNA)
  CNAs = gtools::mixedsort(CNAs)
  nclusters = length(unique(clustering_assignments$cluster))
  ret <- lapply(CNAs, function(s) plot_fit_density_aux_mod(data_hist, 
                                                           densities, s, x, position, nclusters, palette_name))
  return(ret)
}


plot_fit_density_aux_mod = function (df, densities, segment, x, position, nclusters, palette_name) 
{
  df <- df %>% filter(segment_id == segment)

    if (nclusters > 9) {
      getPalette = grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))
      cols = getPalette(nclusters)
  } else {
    cols <- RColorBrewer::brewer.pal(nclusters, palette_name)
  }


  densities <- densities %>% filter(segment_id == segment)
  # cols <- RColorBrewer::brewer.pal(9, palette_name)
  CNA <- df %>% arrange(cluster) %>% dplyr::select(cluster, CNA) %>% 
    unique() %>% pull(CNA)
  clts <- sort(df$cluster %>% unique())
  mu <- plyr::ddply(df, c("modality", "cluster"), summarise, grp.mean=mean(value))
  
  p1 <- ggplot() + 
    geom_histogram(aes(x = value, 
                       fill = factor(cluster, levels = clts)), 
                   bins = 50, 
                   data = df, 
                   position=position,
                   #color = "black", 
                   alpha = 0.6) + 
    geom_vline(data=mu, aes(xintercept=grp.mean, color=factor(cluster, levels = clts)),
               linetype="dashed", alpha=1, size=1)+
    scale_color_manual("Cluster", values = cols, drop = FALSE) +
    facet_wrap(segment_id ~ modality, scales = "free") + 
    guides() + 
    theme_linedraw(base_size = 9) + 
    scale_fill_manual("CN value", values = cols, labels = CNA, drop = FALSE) + 
    labs(x = "Input", y = "Observations") + 
    theme(strip.text.y.right = element_text(angle = 0), 
          axis.text.y = element_blank(), 
          axis.ticks.y = element_blank(), 
          axis.ticks = element_blank())
  p2 <- ggplot() + 
    geom_line(aes(x = X, 
                  y = value, 
                  color = factor(cluster, levels = clts)), 
              data = densities, 
              size = 0.8) + 
    geom_vline(data=mu, aes(xintercept=grp.mean, color=factor(cluster, levels = clts)),
               linetype="dashed", alpha=1, size=1)+
    
    #scale_color_manual("Cluster", values = cols, drop = FALSE) +
    geom_point(aes(x = X, 
                   y = value, 
                   color = factor(cluster, levels = clts)), 
               data = densities, size = 0.3) + 
    facet_wrap(segment_id ~ modality, scales = "free") + 
    geom_text(data=mu, mapping=aes(x=grp.mean, y=0,
                                   label=format(grp.mean, trim=T, digits=3)), 
              size=4, angle=90, vjust=-0.4, hjust=0) +
    labs(title = x$description, subtitle = "Cell cluster assignments and densities") + 
    guides(fill = FALSE) + theme_linedraw(base_size = 9) + 
    scale_color_manual("Clusters", values = cols, drop = FALSE) + 
    labs(x = "Input", y = "Density") + 
    theme(strip.text.y.right = element_text(angle = 0), 
          axis.text.y = element_blank(), 
          axis.ticks.y = element_blank())
  
  cowplot::plot_grid(p2, p1, align = "hv", axis = "lrtb", 
                     ncol = 1)
}


plot_data_mod = function(x,
                         what = 'histogram',
                         ...)
{
  x %>% Rcongas:::sanitize_obj()
  
  if (what == 'histogram')
    return(x %>% plot_data_histogram_mod(...))
  
  if (what == 'lineplot')
    return(x %>% plot_data_lineplot(...))
  
  if (what == 'heatmap'){
    # This returns a more complex assembly
    all_plots = x %>% plot_data_heatmap(...)
    return(all_plots$figure)
  }
  
  if (what == 'mapping')
    return(x %>% plot_data_mapping())
  
  stop("Unrecognised 'what': use any of 'histogram', 'lineplot', 'heatmap' or 'mapping'.")
}

# This function takes in input a congas+ object with a x$input$metadata field, which
# is a dataframe with additional information that may be used for plotting
# In case a user wishes to plot only counts or clusters, the object is not required to have
# the metadata field.
# If to_plot = 'clusters', then this function returns a list of plots, one per segment.
# Each plot is a histogram, colored based on the cluster. The legend indicates both 
# the cluster name and the copy number associated to it.
# single_segment_mode is a boolean. If TRUE, a list of plots (one for each segment) will be returned.
# Otherwise, one plot faceted based on segments and modality is returned. Note that in case
# to_plot = 'clusters' this has no effect on the final plot.
plot_data_histogram_mod = function(x,
                                   to_plot = NULL, 
                                   segments = get_input(x, what = 'segmentation') %>% pull(segment_id),
                                   whichfacet = ggplot2::facet_wrap, bins = 60, position = 'stack',
                                   highlights = FALSE,
                                   single_segment_mode = FALSE, palette_name = NULL, colors = NULL)
{
  if (is.null(to_plot)) {
    to_plot = 'modality'
  } else {
    if (! to_plot %in% colnames(x$input$metadata) & to_plot != 'clusters') {
      stop('Error, to_plot must be either the name of a metadata column or "clusters". Exiting')
    }
  }
  stats_data = Rcongas:::stat(x, what = 'data')
  
  # Data stats
  subtitle = paste0(
    "RNA (",
    stats_data$ncells_RNA,
    " cells, ",
    stats_data$rna_genes,
    " genes) ATAC (",
    stats_data$ncells_ATAC,
    " cells, ",
    stats_data$atac_peaks,
    " peaks)"
  )
  
  # What to print, enhanced facets
  what_rna_lik = case_when(
    stats_data$rna_dtype == "NB" ~ "Negative Binomial",
    stats_data$rna_dtype == "P" ~ "Poisson",
    stats_data$rna_dtype == "G" ~ "Gaussian"
  )
  
  what_atac_lik = case_when(
    stats_data$atac_dtype == "NB" ~ "Negative Binomial",
    stats_data$atac_dtype == "P" ~ "Poisson",
    stats_data$atac_dtype == "G" ~ "Gaussian"
  )
  
  # Normalisation factors
  norm_factors = get_input(x, what = 'normalisation')
  
  # RNA_data values
  what_RNA = get_input(x, what = 'data') %>%
    filter(modality == "RNA")
  
  if (nrow(what_RNA) > 0  && stats_data$rna_dtype != "G")
    what_RNA = Rcongas:::normalise_modality(what_RNA, norm_factors %>% filter(modality == "RNA"))
  
  # ATAC_data values
  what_ATAC = get_input(x, what = 'data') %>%
    filter(modality == "ATAC")
  
  if (nrow(what_ATAC) > 0 && stats_data$atac_dtype != "G")
    what_ATAC = Rcongas:::normalise_modality(what_ATAC, norm_factors %>% filter(modality == "ATAC"))
  
  # Rescale against average normalization factors
  what_ATAC = what_ATAC %>%  mutate(value = value * mean(norm_factors %>% filter(modality == "ATAC") %>%  pull(normalisation_factor)))
  
  what_RNA = what_RNA %>%  mutate(value = value * mean(norm_factors %>% filter(modality == "RNA") %>%  pull(normalisation_factor)))
  
  what = bind_rows(what_RNA, what_ATAC) %>%
    filter(segment_id %in% segments) %>%
    mutate(modality = case_when(
      modality == "RNA" ~ paste(modality, '(', what_rna_lik, ')'),
      modality == "ATAC" ~ paste(modality, '(', what_atac_lik, ')')
    ))
  
  what$modality <- sapply(what$modality %>% strsplit(" "), 
                          function(y) y[1])
  
  # Set some properties of the plot
  if (position != 'stack') alpha = 0.8 else alpha = 1
  
  # Get segments to plot
  if (highlights) {
    CNAs = get_fit(x, "CNA")
    nclusters = CNAs$cluster %>% unique %>% length()
    CNAs = CNAs %>% group_by(segment_id, value) %>% mutate(grp_size = n()) %>% 
      filter(grp_size != nclusters) %>% pull(segment_id) %>% 
      unique
    cli::cli_alert("Plotting segments where different CNAs are present: {.field {CNAs}}.")
  } else {
    CNAs = x$input$segmentation %>% pull(segment_id) %>% unique
    cli::cli_alert("Showing all segments (this plot can be large).")
  }
  
  # Possible behaviours of this function: 
  # 1. to_plot == 'clusters'. In this case you need to get the cluser assignments from the object and then you return a list with one plot per element, colored according to cluster assignments.
  if (to_plot == 'clusters') {
    what = dplyr::left_join(what, get_fit(x, what = "cluster_assignments") %>% dplyr::select(-modality))

    
    if (length(CNAs) == 0) {
      return(ggplot() + geom_blank())
    }
    
    CNA <- Rcongas:::get_CNA(x)
    colnames(CNA)[3] <- "CNA"
    what <- dplyr::left_join(what, CNA)
    what$clusterCNA = paste0(what$cluster, " (CN = ", what$CNA, ")")
    
    if (!is.null(colors)) {
      cols = colors
    } else {
      palette_name = if (!is.null(palette_name)) palette_name else 'Set1'
      nclusters = CNA$cluster %>% unique %>% length()
        if (nclusters > 9) {
          getPalette = grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, palette_name))
          cols = getPalette(nclusters)
        } else {
          cols <- RColorBrewer::brewer.pal(nclusters, palette_name)
        }
    }



    CNAs = gtools::mixedsort(CNAs)
    ret <- sapply(CNAs,  simplify = F, USE.NAMES = T, function(s) {
      ggplot() + 
        geom_histogram(aes_string("value", fill = 'clusterCNA'), 
                       bins = bins, 
                       data = what %>% filter(segment_id == s), 
                       position=position,
                       alpha = alpha) +
                       #color = 'black',
                       #linewidth = 0.1) + 
        #scale_color_manual("Cluster", values = cols, drop = FALSE) +
        facet_wrap(segment_id ~ modality, scales = "free") + 
        guides() + 
        theme_light(base_size = 9) + 
        theme(strip.text = element_text(colour = 'black')) +
        scale_fill_manual("Cluster ", values = cols) + 
        labs(x = "Input", y = "Observations") + 
        theme(strip.text.y.right = element_text(angle = 0), 
              # axis.text.y = element_blank(), 
              # axis.ticks.y = element_blank(), 
              # axis.ticks = element_blank(), 
              legend.position="top",
              legend.direction='vertical') 
    })
    return(ret)
  } 
  
  # 2. to_plot == any metadata column. IN this case you need to get the values to plot from the metadata column. 
  # 3. to_plot == 'modality'. In this case, we don't need to extract anything from the metadata, the variable what 
  # already contains the data we need for plotting.
  
  # Here we are putting inside the variable what the data we need for plotting (so we are finishing preparing 'what')
  # In case to_plot is not null we take its value from the metadata dataframe
  # In case to_plot is 'modality', then there is no need for the object to
  # have a metadata field.
  if (to_plot != 'modality') {
    what = dplyr::left_join(what, x$input$metadata %>% dplyr::select(cell, to_plot))
    what[[to_plot]] = factor(what[[to_plot]], levels = gtools::mixedsort(what[[to_plot]] %>% unique))
  } 
  
  what = what %>%
    mutate(segment_id = factor(segment_id, levels = gtools::mixedsort(segment_id %>% unique))) 
  
  # Two possible behaviours: if single_segment_mode == TRUE, return only one plot faceted. Otherwise return a list of plots.
  
  plot_list = list()
  if (single_segment_mode) {
    plot_list = sapply(CNAs, simplify = F, USE.NAMES = T, function(s) plot_data_histogram_mod_aux(what, s, to_plot, alpha, position, palette_name, colors))
    # for (seg in CNAs) {
    
    #   plot_list[[seg]] = ggplot
    #   (what %>% filter(segment_id == seg))
    # }
  } else {
    plot_list[[1]] = plot_data_histogram_mod_aux(what, s = NULL, to_plot = to_plot, alpha = alpha, position = position, palette_name, colors) + 
      labs(title = x$description, subtitle = subtitle)
  }
  
  # for (pl in seq(1,length(plot_list))) {
  #   p = plot_list[[pl]] +
  #     geom_histogram(aes_string("value", fill = to_plot), alpha=alpha, bins = 50, position=position) +
  #     facet_wrap(segment_id~modality, scales = 'free') +
  #     # facet_wrap(chr ~ modality, scales = 'free') +
  #     theme_light(base_size = 9) + 
  #     labs(x = "Input",
  #          y = 'Observations') +
  #     theme(strip.text.y.right = element_text(angle = 0), 
  #           legend.position="top", legend.direction='horizontal',
  #           legend.title = element_blank()) +
  #     theme(strip.text = element_text(colour = 'black')) 
  
  #   if (to_plot == 'modality') {
  #     p = p + scale_fill_manual(values = cols) 
  #   }
  #   plot_list[[pl]] = p 
  # }
  
  if (length(plot_list) == 1) {
    return(plot_list[[1]])
  } else 
    return(plot_list)
}

plot_data_histogram_mod_aux = function(what, s, to_plot, alpha, position, palette_name, colors = NULL) {
  
  if (!is.null(s)) {
    what = what %>% filter(segment_id == s)
  }
  
  p = ggplot(what) + 
    geom_histogram(aes_string("value", fill = to_plot), 
      alpha=alpha, 
      bins = 50, 
      position=position) +
      #color = 'black', linewidth = 0.1) +
    facet_wrap(segment_id ~ modality, scales = 'free') +
    theme_light(base_size = 9) + 
    labs(x = "Input",
         y = 'Observations') +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 4, min.n = 2)) +
		scale_y_continuous(breaks = scales::pretty_breaks(n = 4, min.n = 2)) +
    theme(strip.text.y.right = element_text(angle = 0), 
          legend.position="top", legend.direction='horizontal',
          legend.title = element_blank()) +
    theme(strip.text = element_text(colour = 'black')) 

  if (!is.null(palette_name)) {
    cols <- RColorBrewer::brewer.pal(length(unique(what[,to_plot])), palette_name)
    p = p + scale_fill_manual(values = cols) 
  } else if(!is.null(colors)) {
    p = p + scale_fill_manual(values = colors) 
  }
  
  if (to_plot == 'modality') {
    cols = Rcongas:::modality_colors(what$modality %>% unique)
    p = p + scale_fill_manual(values = cols) 
  } 
  return(p)
}

