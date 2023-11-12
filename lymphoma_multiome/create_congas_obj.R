
library(Rcongas)
library(stringi)
library(dplyr)
library(ggplot2)
library(stringr)

script_folder = './'
data_folder = 'data/'

source(paste0(script_folder, "/utils/congas_plot_mod.R"))
source(paste0(script_folder, "/utils/outliers_new_function.R"))
source(paste0(script_folder,'/utils/create_tibble.R'))



setwd(paste0(data_folder, "/10x_multiome"))

sample = 'lymphoma'
lik_rna = 'NB'
lik_atac = 'NB'
prefix = paste0(sample, "/bimodal_rna_", lik_rna, '_atac_', lik_atac, "_new/")

out.dir = paste0(prefix, "/congas_data/")
fig.dir = paste0(prefix, "/congas_figures/")

if (!dir.exists(out.dir)) {dir.create(out.dir, recursive = T)}
if (!dir.exists(fig.dir)) {dir.create(fig.dir, recursive = T)}
# First create a dataset with only normal cells (T and non-Tcells) 
atac = readr::read_tsv(paste0(sample, '/atac/counts_final.tsv'),
                       col_types = "cdcii",
                       col_names = c('cell', 'value', 'chr', 'from', 'to')) 
atac = atac[-c(1),] 
atac = atac %>% mutate(cell = paste0(cell, '-ATAC'))
norm_atac = Rcongas:::auto_normalisation_factor(atac) %>%
  mutate(modality = 'ATAC')

rna = readr::read_tsv(paste0(sample, '/rna/counts_final.tsv'),
                      col_types = "ciccii",
                      col_names = c("cell", "value",  'gene', 'chr', 'from', 'to'))
rna = rna[-c(1),] 
rna = rna %>% mutate(cell = paste0(cell, '-RNA'))

rna$gene[is.na(rna$gene)] = '' 
norm_rna = Rcongas:::auto_normalisation_factor(rna) %>%
  mutate(modality = 'RNA')


rna = rna %>% filter_known_genes(what='r')
# 
all_genes = rna$gene %>% unique
mito = all_genes %>% str_starts(pattern = 'MT-')
all_genes = setdiff(all_genes, all_genes[mito])
rna = rna %>% dplyr::filter(gene %in% all_genes)


segments = readr::read_delim(paste0(data_folder, 'chromosome_arms_hg38.tsv'),
                             col_types = "cii", delim = ' ')  %>% mutate(copies = as.integer(2))




# Remove these chromosomes
segments = segments %>% dplyr::filter(chr != 'chrX', chr != 'chrY')
atac = atac %>% mutate(value = as.integer(value))
rna = rna %>% mutate(value = as.integer(value))


# 
# DATASET 1 
# 
library(tidyverse)


x = init(
  rna = rna,
  atac = atac,
  segmentation = segments, #%>% mutate(copies = as.integer(round(copies))),
  rna_normalisation_factors = norm_rna,
  atac_normalisation_factors = norm_atac,
  rna_likelihood = lik_rna, 
  atac_likelihood = lik_atac,
  description = paste0('Bimodal ', sample))



# Now plot the number of nonzero cells and the number of peaks per segment
ggplot(x$input$segmentation, aes(x=ATAC_nonzerocells, y=ATAC_peaks)) + geom_point()
ggplot(x$input$segmentation, aes(x=RNA_nonzerocells, y=RNA_genes)) + geom_point()

sortedATAC = x$input$segmentation %>% arrange(ATAC_nonzerocells)
sortedRNA = x$input$segmentation %>% arrange(RNA_nonzerocells)

segments = x$input$segmentation 


ATAC_threshold =12000
RNA_threshold = 12000
segmentsRemove = segments %>% filter((RNA_nonzerocells < !!RNA_threshold) | (ATAC_nonzerocells < !!ATAC_threshold))
segmentsKeep = segments %>% filter((RNA_nonzerocells >= !!RNA_threshold) & (ATAC_nonzerocells >= !!ATAC_threshold))

x = Rcongas:::select_segments(x, segment_ids = segmentsKeep$segment_id)


# Filter post mapping segments data
s_q = c(0.01, .99)
x = x %>% Rcongas:::filter_outliers(lower_quantile = s_q[1], upper_quantile = s_q[2], frequency_cutoff = 0)# , frequency_cutoff = 0)#, action = 'cap')


ggsave(
  paste0(fig.dir, "3.1.Histogram_longest_20_segments.png"),
  plot_data(
    x,
    what = 'histogram',
    segments = get_input(x, what = 'segmentation') %>%
      mutate(L = to - from) %>%
      dplyr::arrange(dplyr::desc(L)) %>%
      top_n(20) %>%
      pull(segment_id)
  ),
  width = 18,
  height = 18
)
ggsave(
  paste0(fig.dir, "3.1.Histogram_all_segments_sameOut.png"),
  plot_data(
    x,
    what = 'histogram',
    segments = get_input(x, what = 'segmentation') %>%
      pull(segment_id)
  ),
  width = 18,
  height = 18
)

