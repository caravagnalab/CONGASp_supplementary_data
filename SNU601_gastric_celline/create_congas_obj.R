
library(Rcongas)
library(stringi)
library(dplyr)
library(ggplot2)
library(stringr)

script_folder = './'
data_folder = 'data/'

source(paste0(script_folder, "/utils/congas_plot_mod.R"))
source(paste0(script_folder, "/utils/outliers_new_function.R"))

setwd(paste0(data_folder, "/SNU601/"))

samples = list.files('atac', pattern = 'counts_final.tsv')

lik_rna = 'NB'
lik_atac = 'NB'
prefix = paste0("bimodal_rna_", lik_rna, '_atac_', lik_atac, "/")

out.dir = paste0(prefix, "/congas_data/")
fig.dir = paste0(prefix, "/congas_figures/")

if (!dir.exists(out.dir)) {dir.create(out.dir, recursive = T)}
if (!dir.exists(fig.dir)) {dir.create(fig.dir, recursive = T)}



atac = readr::read_tsv(paste0('atac/', samples[1]),
                       col_types = "cdcii",
                       col_names = c('cell', 'value', 'chr', 'from', 'to')) 
atac = atac[-c(1),] 
atac = atac %>% mutate(cell = paste0(cell, '-atac'))

samples_rna = list.files('rna', pattern = 'counts_final.tsv')
rna = readr::read_tsv(paste0('rna/', samples_rna[1]),
                      col_types = "ciccii",
                      col_names = c("cell", "value", 'gene', 'chr', 'from', 'to'))
rna = rna[-c(1),] 
rna = rna %>% mutate(cell = paste0(cell, '-', 'rna'))

norm_atac = Rcongas:::auto_normalisation_factor(atac) %>%
  mutate(modality = 'ATAC')

norm_rna = Rcongas:::auto_normalisation_factor(rna) %>%
  mutate(modality = 'RNA')

segments = readRDS('atac/segmentation.rds')

# res_rna = Rcongas:::clean_outlers_persegment(segments, 'RNA', rna, norm_rna)
atac = atac %>% mutate(value = as.integer(value))
library(tidyverse)

norm_atac = Rcongas:::auto_normalisation_factor(atac) %>%
  mutate(modality = 'ATAC')

rna = rna %>% filter_known_genes(what='r')
#
all_genes = rna$gene %>% unique
mito = all_genes %>% str_starts(pattern = 'MT-')
all_genes = setdiff(all_genes, all_genes[mito])
rna = rna %>% dplyr::filter(gene %in% all_genes)
# Remove these chromosomes
segments = segments %>% dplyr::filter(chr != 'chrX', chr != 'chrY')



x = init(
  rna = rna  %>% select(chr, from, to, cell, value, gene) %>% 
    replace_na(replace = list(gene= '', chr = '', from = 0, to = 0)),
  atac =  atac %>% select(chr, from, to, cell, value),
  segmentation = segments %>% mutate(copies = as.integer(copies), from = as.integer(from), to = as.integer(to)),
  rna_normalisation_factors = norm_rna,
  atac_normalisation_factors = norm_atac,
  rna_likelihood = lik_rna, 
  atac_likelihood = lik_atac,
  description = paste0('Bimodal SNU-601'))



sortedATAC = x$input$segmentation %>% arrange(ATAC_nonzerocells)
sortedRNA = x$input$segmentation %>% arrange(RNA_nonzerocells)



segments = x$input$segmentation %>% filter(chr != 'chr13') 


x = Rcongas:::select_segments(x, segment_ids = segments$segment_id)


# Filter post mapping segments data
s_q = c(0.01, .99)
x = x %>% Rcongas:::filter_outliers(lower_quantile = s_q[1], upper_quantile = s_q[2], frequency_cutoff = 0)# , frequency_cutoff = 0)#, action = 'cap')

#########################
ggsave(
  paste0(fig.dir, "3.1.Histogram_longest_20_segments_0outN.png"),
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
  paste0(fig.dir, "3.1.Histogram_all_segments_0out.png"),
  plot_data(
    x,
    what = 'histogram',
    segments = get_input(x, what = 'segmentation') %>%
      pull(segment_id)
  ),
  width = 18,
  height = 18
)


saveRDS(x, paste0(out.dir, "5.rcongas_noOutliers.rds"))



