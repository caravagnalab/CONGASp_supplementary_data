library(Rcongas)
library(stringi)
library(dplyr)
library(ggplot2)
library(stringr)
library(tidyverse)

script_folder = './'
data_folder = 'data/'


script_folder = './'
data_folder = 'data/'

source(paste0(script_folder, "/utils/congas_plot_mod.R"))
source(paste0(script_folder, "/utils/outliers_new_function.R"))

setwd(paste0(data_folder, "/celline_prostate/"))

samples = list.files('atac', pattern = 'counts_final.tsv')

lik_rna = 'NB'
lik_atac = 'NB'
prefix = paste0('all', "/bimodal_rna_", lik_rna, '_atac_', lik_atac, "_new/")

out.dir = paste0(prefix, "/congas_data/")
fig.dir = paste0(prefix, "/congas_figures/")

if (!dir.exists(out.dir)) {dir.create(out.dir, recursive = T)}
if (!dir.exists(fig.dir)) {dir.create(fig.dir, recursive = T)}

# First create a dataset with only normal cells (T and non-Tcells) 
atac = tibble(cell = character(), value = numeric(), chr = character(), from = integer(), to = integer())
metadata_atac = tibble::tibble(cell = character(), type = character())
for (i in seq(1, length(samples))){  
  sample = strsplit(samples[i], '_')[[1]][2] 
  atac_current = readr::read_tsv(paste0('atac/', samples[i]),
                         col_types = "cdcii",
                         col_names = c('cell', 'value', 'chr', 'from', 'to')) 
  atac_current = atac_current[-c(1),] 
  atac_current = atac_current %>% mutate(cell = paste0(cell, !!sample))
  atac = bind_rows(atac, atac_current)
  metadata_current = atac_current %>% select(cell) %>% unique %>% mutate(type = !!sample )
  metadata_atac = bind_rows(metadata_atac, metadata_current)
}
samples_rna = list.files('rna', pattern = 'counts_final.tsv')
rna = tibble(gene = character(),  chr = character(), from = integer(), to = integer(), cell = character(), value = integer())
metadata_rna = tibble::tibble(cell = character(), type = character())
for (i in seq(1, length(samples))){ 
  sample = strsplit(samples_rna[i], '_')[[1]][2] 
  rna_current = readr::read_tsv(paste0('rna/', samples_rna[i]),
                                col_types = "ciccii",
                                col_names = c("cell", "value", 'gene', 'chr', 'from', 'to'))
  rna_current = rna_current[-c(1),] 
  rna_current = rna_current %>% mutate(cell = paste0(cell, '-', !!sample))
  rna = bind_rows(rna, rna_current)
  metadata_current = rna_current %>% select(cell) %>% unique %>% mutate(type = !!sample)
  metadata_rna = bind_rows(metadata_rna, metadata_current)
}
norm_atac = Rcongas:::auto_normalisation_factor(atac) %>%
  mutate(modality = 'ATAC')

norm_rna = Rcongas:::auto_normalisation_factor(rna) %>%
  mutate(modality = 'RNA')


segments_hg19 = readr::read_delim(paste0('cellline_segments.tsv'),
                                  col_types = "cii", delim = '\t')

lifted_start = readr::read_delim(paste0('lifted_start.bed'),
                                 col_types = "cii", delim = '\t', 
                                 col_names = c('chr', 'from', 'to', 'old', 'n'))

lifted_end = readr::read_delim(paste0('lifted_end.bed'),
                                 col_types = "cii", delim = '\t', 
                                 col_names = c('chr', 'from', 'to', 'old', 'n'))

segments = readr::read_delim(paste0('cellline_segments.tsv'),
                             col_types = "cii", delim = '\t') 
segments$from = lifted_start$from
segments$to = lifted_end$from
segments = segments %>% select(chr, from, to, copies)

res_rna = Rcongas:::clean_outlers_persegment(segments, 'RNA', rna, norm_rna)
atac = atac %>% mutate(value = as.integer(value))
library(tidyverse)
res_atac = Rcongas:::clean_outlers_persegment(segmentation = segments, modality = 'ATAC', 
                                              data = atac, normalisation_factors = norm_atac)

norm_atac = Rcongas:::auto_normalisation_factor(atac) %>%
  mutate(modality = 'ATAC')

norm_rna = Rcongas:::auto_normalisation_factor(rna) %>%
  mutate(modality = 'RNA')

rna = res_rna$data_cleaned %>% filter(!is.na(chr))
atac = res_atac$data_cleaned %>% replace(is.na(.), '')

rna = rna %>% filter_known_genes(what='r')
#
all_genes = rna$gene %>% unique
mito = all_genes %>% str_starts(pattern = 'MT-')
all_genes = setdiff(all_genes, all_genes[mito])
rna = rna %>% dplyr::filter(gene %in% all_genes)

# Remove these chromosomes
segments = segments %>% dplyr::filter(chr != 'chrX', chr != 'chrY')



x = init(
  rna = rna %>% select(chr, from, to, cell, value, gene),
  atac = atac %>% select(chr, from, to, cell, value),
  segmentation = segments %>% mutate(copies = as.integer(round(copies))),
  rna_normalisation_factors = norm_rna,
  atac_normalisation_factors = norm_atac,
  rna_likelihood = lik_rna, 
  atac_likelihood = lik_atac,
  description = paste0('Bimodal ', sample))



ggsave(
  paste0(fig.dir, "1.1.Histogram_longest_20_segments.png"),
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
  paste0(fig.dir, "1.1.Histogram_all_segments.png"),
  plot_data(
    x,
    what = 'histogram',
    segments = get_input(x, what = 'segmentation') %>%
      pull(segment_id)
  ),
  width = 18,
  height = 18
)

nggsave(
  paste0(fig.dir, "1.2.Cell_heatmap.png"),
  plot_data(
    x,
    what = 'heatmap',
    segments = get_input(x, what = 'segmentation') %>%
      mutate(L = to - from) %>%
      arrange(dplyr::desc(L)) %>%
      #top_n(20) %>%
      pull(segment_id)
  ),
  width = 12,
  height = 9
)
x$input$metadata = rbind(metadata_atac, metadata_rna) %>% filter(cell %in% unique(x$input$dataset$cell))
# Modificato plot_data in mdodo da plottare gli istogrammi colorati per celltype
ggsave(paste0(fig.dir, "1.1.Histogram_type.png"),
       plot_data_mod(x, to_plot = 'type', position = 'stack'),
       width = 19, height = 19)
ggsave(
  paste0(fig.dir, "1.3.Events_mapping.png"),
  plot_data(x, what = 'mapping'),
  width = 4,
  height = 16
)

sortedATAC = x$input$segmentation %>% arrange(ATAC_nonzerocells)
sortedRNA = x$input$segmentation %>% arrange(RNA_nonzerocells)



segments = x$input$segmentation 

segmentsRemove = segments %>% filter(RNA_genes < 100)
segmentsKeep = segments %>% filter(RNA_genes >= 100

xdiscarded = Rcongas:::select_segments(x, segment_ids = segmentsRemove$segment_id)

ggsave(
  paste0(fig.dir, "1.1.Histogram_discarded_segments.png"),
  plot_data(
    xdiscarded,
    what = 'histogram',
    segments = get_input(xdiscarded, what = 'segmentation') %>%
      pull(segment_id)
  ),
  width = 18,
  height = 18
)

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
  paste0(fig.dir, "3.1.Histogram_all_segments.png"),
  plot_data(
    x,
    what = 'histogram',
    segments = get_input(x, what = 'segmentation') %>%
      pull(segment_id)
  ),
  width = 18,
  height = 18
)

resistance = data.frame(type = unique(x$input$metadata$type), resistance = c('sensitive', 'sensitive', 'resistant', 'resistant'))

x$input$metadata = x$input$metadata %>% filter(cell %in% unique(x$input$dataset$cell)) %>%
  left_join(resistance)


# Modificato plot_data in mdodo da plottare gli istogrammi colorati per celltype
ggsave(paste0(fig.dir, "3.1.Histogram_type_0out.png"),
       plot_data_mod(x, to_plot = 'type', position = 'stack'),
       width = 19, height = 19)

ggsave(paste0(fig.dir, "3.1.Histogram_resistance.png"),
       plot_data_mod(x, to_plot = 'resistance', position = 'stack'),
       width = 19, height = 19)

saveRDS(x, paste0(out.dir, "5.rcongas_noOutliers.rds"))

