#'@parameter action. Can be one of c('median', 'minMax') 
cap_outliers = function(x, lower_quantile = 0.01, upper_quantile = 0.99, action = 'minMax'){ 
  retained_rna = retained_atac = norm_factors_rna = norm_factors_atac = NULL
  
  if(Rcongas:::has_rna(x))
  {
    retained_rna = Rcongas:::compute_outliers_per_segment(x, modality = "RNA", lower_quantile, upper_quantile) 
    # %>% mutate(outlier = as.character(outlier)) %>%
    #   mutate(outlier = case_when(value < q_min ~ 'LEFT',
    #                              value > q_max ~ 'RIGHT',
    #                              TRUE ~ 'FALSE'))
    median_cells_signal = retained_rna %>% filter(!outlier) %>% group_by(segment_id) %>% summarise(median_value = median(value), max_value= max(value), min_value = min(value))
    
    if (x$input$dataset$value_type[1] == 'G' ){ 
       retained_rna = retained_rna %>% left_join(median_cells_signal) %>% 
        mutate(value = ifelse(outlier, median_value, value)) %>%
        # case_when(outlier == 'RIGHT' ~ max_value,
        #                 outlier == 'LEFT' ~ min_value,
        #                 TRUE ~ value)) %>%
        # 
        dplyr::select(segment_id, cell, value, modality, value_type)
    } else{ 
      retained_rna = retained_rna %>% left_join(Rcongas:::get_input(x, what = 'normalisation'))
      retained_rna = retained_rna %>% left_join(median_cells_signal) %>% 
        mutate(original_value = ifelse(outlier, round(median_value * normalisation_factor), original_value)) %>%
        dplyr::select(segment_id, cell, original_value, modality, value_type) %>% rename(value = original_value)
    }
  }
  
  if(Rcongas:::has_atac(x))
  {
    retained_atac = Rcongas:::compute_outliers_per_segment(x, modality = "ATAC", lower_quantile, upper_quantile) 
    # %>% mutate(outlier = as.character(outlier)) %>%
    #   mutate(outlier = case_when(value < q_min ~ 'LEFT',
    #                              value > q_max ~ 'RIGHT',
    #                              TRUE ~ 'FALSE'))
    retained_atac = retained_atac %>% left_join(Rcongas:::get_input(x, what = 'normalisation'))
    # In any segment, if x is an outlier I want to replace it with x/lib_size = max(y/libsizey) Where y are all cells not outliers
    #So the normalized value of x will be x = max(y/libsizey)*libsize(x). I already have normalized ys
    median_cells_signal = retained_atac %>% filter(!outlier) %>% group_by(segment_id) %>% summarise(median_value = median(value), max_value = max(value), min_value = min(value))
    retained_atac = retained_atac %>% left_join(median_cells_signal) %>% 
      mutate(original_value = ifelse(outlier, round(median_value * normalisation_factor), original_value)) %>%
      # case_when(outlier == 'LEFT' ~ round(min_value * normalisation_factor),
      #                          outlier == 'RIGHT' ~ round(max_value * normalisation_factor),
      #                          TRUE ~ original_value)) %>%
      # 
      dplyr::select(segment_id, cell, original_value, modality, value_type) %>% rename(value = original_value)
    # norm_factors_atac = aux_fun_nf(data = retained_atac, modality = "ATAC")
  }
  
  retained = bind_rows(retained_rna, retained_atac)
  # norm_factors = bind_rows(norm_factors_rna, norm_factors_atac)
  
  # Copy input x object, before making a new one
  # x_original = x
  
  # Rebuild loses some info, not super important
  x$input$dataset = retained
  # x$input$normalisation = norm_factors
  
  # x$input$normalisation = x %>% 
  #   get_input(what = 'normalisation') %>% 
  #   filter(cell %in% retained$cell %>% unique())
  
  x$description = paste0(
    x$description, 
    "; post-map q[", lower_quantile, ', ', upper_quantile, '] ', 'cap')
  
  return(x)
  
}


# replace_outliers = function(retained_data, modality, )

