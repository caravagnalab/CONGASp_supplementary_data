#' This function takes in input a SummarizedExperiment object with rowranges and  and creates a tibble with counts needed to create the Rcongas object.
#' @param counts (Required). Count matrix with cells in columns and genes/bins in rows. 
#' @param modality. It has to be either ATAC or RNA based on the content of the count matrix.
#' @param output_dir. Directory where to save the tibble. 
#' @param features Dataframe containing columns chr, from, to and (only if RNA) gene. 
#' @param file_name. Name to give to the counts file. Default is counts_final.tsv  
create_congas_tibble = function(counts, modality = 'ATAC', save_dir = getwd(),
                                features = NULL,
                                output_file = 'counts_final.tsv') {
  # counts = as.matrix(assay(sumExp))
  # counts  = SummarizedExperiment::assay(sumExp)
  cli::cli_alert_info("Saving counts\n")
  writeMM(as(counts, "sparseMatrix"), paste0(save_dir, "/counts.mtx"))
  
  # Now create tibble cell chr from to value  
  counts = readr::read_delim(paste0(save_dir, "/counts.mtx"),
                             comment = '%', delim = ' ') 
  
  colnames(counts) =  c('feature_index', 'cell_index', 'value')
  
  if (!is.null(features)){
    features = features %>% mutate(feature_index = seq(1, nrow(.)))
  } 
  
  if (modality == 'ATAC' & is.null(features)) {
    stop('Please provide a dataframe with c("chr", "from", "to") through the parameter features')
    # features = data.frame(feature = rownames(sumExp)) %>% mutate(feature_index = seq(1, nrow(.))) %>%
    #   tidyr::separate(col = 'feature', sep = ':', into = c('chr', 'coordinates')) %>%
    #   tidyr::separate(col = 'coordinates', sep ='-', into = c('from', 'to'))
  } else if (modality == 'RNA' & is.null(features)) { 
    stop('Please provide a dataframe with c("chr", "from", "to", "gene") through the parameter features')
    # # For RNA access rowranges in the summarizedexperiment sim object
    # features = as.data.frame(sim@rowRanges) %>% 
    #   mutate(feature_index = seq(1, nrow(.))) %>%
    #   dplyr::select(-c(strand, width)) %>% 
    #   dplyr::rename(chr = seqnames, gene = gene_name, from = start, to = end)
    # rownames(features) = NULL
  }
  cells = data.frame(cell = colnames(counts)) %>% mutate(cell_index = seq(1, nrow(.)))
  
  counts = counts %>% dplyr::left_join(features) %>% dplyr::left_join(cells)
  
  if (modality == 'ATAC') {
    counts = counts %>% dplyr::select(c("cell", "value",  'chr', 'from', 'to')) %>%
      mutate(from = as.integer(from), to = as.integer(to), chr = as.character(chr))
  } else {
    counts = counts %>% dplyr::select(c("cell", "value",  'gene', 'chr', 'from', 'to')) %>%
      mutate(from = as.integer(from), to = as.integer(to))
  }
  
  readr::write_tsv(counts, file = paste0(save_dir, '/', output_file))
  # saveRDS(counts, paste0(save_dir, '/counts_final.rds'))
  
  
}