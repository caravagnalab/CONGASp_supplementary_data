
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

# Apply some basic filters - cap outliers and remove known genes
# atac = filter_counts_by_quantile(atac, upper_quantile = 0.95)
# rna = filter_counts_by_quantile(rna, upper_quantile = 0.95)
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

res = Rcongas:::clean_outlers_persegment(segmentation = segments, modality = 'RNA', 
                                         data = rna, normalisation_factors = norm_rna)

res_atac = Rcongas:::clean_outlers_persegment(segmentation = segments, modality = 'ATAC', 
                                         data = atac, normalisation_factors = norm_atac)
# p = cowplot::plot_grid(plotlist = res$persegment_plot)
# ggsave(paste0(fig.dir, 'segments_ouliers_rna.pdf'), p, width = 15, height = 15)

x = init(
  rna = res$data_cleaned,
  atac = res_atac$data_cleaned,
  segmentation = segments, #%>% mutate(copies = as.integer(round(copies))),
  rna_normalisation_factors = norm_rna,
  atac_normalisation_factors = norm_atac,
  rna_likelihood = lik_rna, 
  atac_likelihood = lik_atac,
  description = paste0('Bimodal ', sample))

# Filter post mapping segments data
# s_q = c(0.01, .99)
# x = x %>% filter_outliers(lower_quantile = s_q[1], upper_quantile = s_q[2])#, action = 'cap')

#x = readRDS("1.rcongas_raw.rds")

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
# x$input$metadata = rbind(metadata_atac, metadata_rna) %>% filter(cell %in% unique(x$input$dataset$cell))
# # Modificato plot_data in mdodo da plottare gli istogrammi colorati per celltype
# ggsave(paste0(fig.dir, "1.1.Histogram_type.png"),
#        plot_data_mod(x, to_plot = 'type', position = 'stack'),
#        width = 19, height = 19)
ggsave(
  paste0(fig.dir, "1.3.Events_mapping.png"),
  plot_data(x, what = 'mapping'),
  width = 4,
  height = 16
)
tmp = x$input$dataset %>% filter(modality == 'ATAC') %>% left_join(x$input$normalisation) %>% 
  mutate(value = (value / normalisation_factor) * mean(x$input$normalisation %>% filter(modality == 'ATAC') %>% pull(normalisation_factor)))
ggplot(tmp %>% filter(segment_id == 'chr7:0:60100000'), aes(x = value)) + geom_histogram(bins = 100)
maybeAmp = tmp %>% filter(segment_id == 'chr7:0:60100000' & value > 300) %>% pull(cell)
maybeAmp = c(maybeAmp, gsub('ATAC', 'RNA', maybeAmp))
x$input$metadata = tibble(cell = x$input$normalisation$cell) %>% mutate(amplified = cell %in% !!maybeAmp)
ggsave(paste0(fig.dir, "1.Histogram_preliminaryAMplification.png"),
       plot_data_mod(x, to_plot = 'amplified', position = 'stack'),
       width = 19, height = 19)

# Rcongas raw
saveRDS(x, paste0(out.dir, "1.rcongas_raw.rds"))
x = readRDS(paste0(out.dir, "1.rcongas_raw.rds"))

# Now plot the number of nonzero cells and the number of peaks per segment
ggplot(x$input$segmentation, aes(x=ATAC_nonzerocells, y=ATAC_peaks)) + geom_point()
ggplot(x$input$segmentation, aes(x=RNA_nonzerocells, y=RNA_genes)) + geom_point()

sortedATAC = x$input$segmentation %>% arrange(ATAC_nonzerocells)
sortedRNA = x$input$segmentation %>% arrange(RNA_nonzerocells)
sortedATAC
sortedRNA

segments = x$input$segmentation 


# segments = segments %>% filter(ATAC_nonzerocells>=4326 & RNA_nonzerocells >=3060)
ATAC_threshold =12000# 7000
RNA_threshold = 12000# 4000
segmentsRemove = segments %>% filter((RNA_nonzerocells < !!RNA_threshold) | (ATAC_nonzerocells < !!ATAC_threshold))
segmentsKeep = segments %>% filter((RNA_nonzerocells >= !!RNA_threshold) & (ATAC_nonzerocells >= !!ATAC_threshold))

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

xdiscarded = xdiscarded %>% cap_outliers(lower_quantile = 0.01, upper_quantile = 0.99)

xdiscarded$input$metadata = x$input$metadata %>% filter(cell %in% unique(xdiscarded$input$dataset$cell))
# Modificato plot_data in mdodo da plottare gli istogrammi colorati per celltype
ggsave(paste0(fig.dir, "1.1.Histogram_discarded_segments_amplified.png"),
       plot_data_mod(xdiscarded, to_plot = 'amplified', position = 'stack'),
       width = 19, height = 19)

x = Rcongas:::select_segments(x, segment_ids = segmentsKeep$segment_id)


# Filter post mapping segments data
s_q = c(0.01, .99)
x_noout = x %>% Rcongas:::filter_outliers(lower_quantile = s_q[1], upper_quantile = s_q[2], frequency_cutoff = 0)# , frequency_cutoff = 0)#, action = 'cap')
x_back = x
x = cap_outliers(x,lower_quantile = s_q[1] , upper_quantile = s_q[2])#  lower_quantile = s_q[1], upper_quantile = s_q[2])

outliers = setdiff(x$input$normalisation$cell, x_noout$input$normalisation$cell)

x$input$metadata = left_join(x$input$metadata, tibble(cell = x$input$normalisation$cell) %>% mutate(is_outlier = cell %in% !!outliers))

ggsave(paste0(fig.dir, "2.1.Histogram_outliers.png"),
       plot_data_mod(x, to_plot = 'is_outlier', position = 'stack'),
       width = 19, height = 19)
# tmp_out = paste0(fig.dir, '/single_segments/')
# dir.create(tmp_out)
# noutliers = list()
# for (seg in x$input$segmentation$segment_id){ 
#   tmp = Rcongas:::select_segments(x, segment_ids = c(seg))
#   tmp_outliers = Rcongas:::filter_outliers(tmp, lower_quantile = 0.01, upper_quantile = 0.99, frequency_cutoff = 0)
#   tmp= cap_outliers(tmp, lower_quantile = 0.01, upper_quantile = 0.99)
#   
#   outliers = setdiff(tmp$input$normalisation$cell, tmp_outliers$input$normalisation$cell)
#   
#   tmp$input$metadata = tibble(cell = tmp$input$normalisation$cell) %>% mutate(is_outlier = cell %in% !!outliers)
#   
#   ggsave(paste0(tmp_out, '/', seg, 'capped.png'),
#          plot_data_mod(tmp, to_plot = 'is_outlier', position = 'stack'),
#          width = 19, height = 19)
#   
#   noutliers[[seg]] = length(outliers) 
# }


ggsave(
  paste0(fig.dir, "2.1.Histogram_longest_20_segments.png"),
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
  paste0(fig.dir, "2.1.Histogram_all_segments.png"),
  plot_data(
    x,
    what = 'histogram',
    segments = get_input(x, what = 'segmentation') %>%
      pull(segment_id)
  ),
  width = 18,
  height = 18
)

ggsave(
  paste0(fig.dir, "2.2.Cell_heatmap.png"),
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

ggsave(
  paste0(fig.dir, "2.3.Events_mapping.png"),
  plot_data(x, what = 'mapping'),
  width = 4,
  height = 16
)
# x$input$metadata = rbind(metadata_atac, metadata_rna) %>% filter(cell %in% unique(x$input$dataset$cell))
# Modificato plot_data in mdodo da plottare gli istogrammi colorati per celltype
ggsave(paste0(fig.dir, "2.1.Histogram_amp.png"),
       plot_data_mod(x, to_plot = 'amplified', position = 'stack'),
       width = 19, height = 19)
saveRDS(x, paste0(out.dir, "2.rcongas_filtered.rds"))


#########################
x_noout_rna = x_noout$input$normalisation %>% filter(modality == 'RNA') %>% pull(cell)
x_noout_atac = x_noout$input$normalisation %>% filter(modality == 'ATAC') %>% pull(cell)

x_noout_rna = gsub('-RNA', '', x_noout_rna)
x_noout_atac = gsub('-ATAC', '', x_noout_atac)
x_noout_cells = intersect(x_noout_atac, x_noout_rna)
x_noout_cells = c(paste0(x_noout_cells, '-ATAC'), paste0(x_noout_cells, '-RNA'))

x_noout$input$dataset = x_noout$input$dataset %>% filter(cell %in% !!x_noout_cells)
x_noout$input$normalisation = x_noout$input$normalisation %>% filter(cell %in% !!x_noout_cells)
ggsave(
  paste0(fig.dir, "3.1.Histogram_longest_20_segments_sameOut.png"),
  plot_data(
    x_noout,
    what = 'histogram',
    segments = get_input(x_noout, what = 'segmentation') %>%
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
    x_noout,
    what = 'histogram',
    segments = get_input(x_noout, what = 'segmentation') %>%
      pull(segment_id)
  ),
  width = 18,
  height = 18
)

x_noout$input$metadata = x$input$metadata %>% filter(cell %in% unique(x_noout$input$dataset$cell))
# Modificato plot_data in mdodo da plottare gli istogrammi colorati per celltype
ggsave(paste0(fig.dir, "3.1.Histogram_amp.png"),
       plot_data_mod(x_noout, to_plot = 'amplified', position = 'stack'),
       width = 19, height = 19)

saveRDS(x_noout, paste0(out.dir, "4.rcongas_noSameOutliers.rds"))



top_segments = get_input(x, what = 'segmentation') %>%
  # mutate(L = to - from) %>%
  dplyr::arrange(dplyr::desc(ATAC_peaks)) %>%
  top_n(20) %>%
  pull(segment_id)

top_x = Rcongas:::select_segments(x_noout, segment_ids = top_segments)


ggsave(
  paste0(fig.dir, "4.1.Histogram_all_segments.png"),
  plot_data(
    top_x,
    what = 'histogram',
    segments = get_input(x_noout, what = 'segmentation') %>%
      pull(segment_id)
  ),
  width = 18,
  height = 18
)

top_x$input$metadata = x$input$metadata %>% filter(cell %in% unique(top_x$input$dataset$cell))
# Modificato plot_data in mdodo da plottare gli istogrammi colorati per celltype
ggsave(paste0(fig.dir, "4.2.Histogram_amplified.png"),
       plot_data_mod(top_x, to_plot = 'amplified', position = 'stack'),
       width = 19, height = 19)
saveRDS(top_x, paste0(out.dir, "3.rcongas_top20.rds"))



chr = 'chr7'
from = 0
to = 60100000

atac_chr7 = atac %>% filter(chr == 'chr7' & from >= !!from & to <= !!to) %>% left_join(norm_atac) %>% mutate(value = (value / normalisation_factor))

atac_chr7 = atac_chr7 %>% left_join(x$input$metadata)

ggplot(atac_chr7)+ geom_histogram(aes(x = value, fill =  amplified), bins = 100)

p = ggplot(atac_chr7 %>% filter(value > 2)) + geom_point(aes(x = from, y = value, color = amplified), size = 0.5)
p
ggsave('chr7.png', plot = p)

