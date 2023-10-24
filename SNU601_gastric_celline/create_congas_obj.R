
library(Rcongas)
library(stringi)
library(dplyr)
library(ggplot2)
library(stringr)

script_folder = './'
data_folder = 'data/'

source(paste0(script_folder, "/utils/congas_plot_mod.R"))
source(paste0(script_folder, "/utils/outliers_new_function.R"))

setwd(paste0(data_folder, "/SNU650/"))

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

# Save rds files for later
saveRDS(rna, paste0(out.dir, "0.rna.rds"))
saveRDS(norm_rna, paste0(out.dir, "0.norm_rna.rds"))
saveRDS(atac, paste0(out.dir, "0.atac.rds"))
saveRDS(norm_atac, paste0(out.dir, "0.norm_atac.rds"))
saveRDS(segments, paste0(out.dir, "0.segments.rds"))
# saveRDS(rbind(metadata_rna, metadata_atac), paste0(out.dir, '0.celltypes.rds'))

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

ggsave(
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

# metadata_atac = tibble(Clone) %>% mutate(cell=names(Clone)) %>% dplyr::rename(type = Clone) 

# metadata_atac$cell = paste0(gsub('-', '.', metadata_atac$cell), 'counts')

x$input$metadata = metadata_atac %>% mutate(type = factor(type))#  rbind(metadata_atac, metadata_rna) %>% filter(cell %in% unique(x$input$dataset$cell))
# Modificato plot_data in mdodo da plottare gli istogrammi colorati per celltype
ggsave(paste0(fig.dir, "1.1.Histogram_type.png"),
       plot_data_mod(x, to_plot = 'type', position = 'stack', palette_name = 'Set1'),
       width = 19, height = 19)
ggsave(
  paste0(fig.dir, "1.3.Events_mapping.png"),
  plot_data(x, what = 'mapping'),
  width = 4,
  height = 16
)

# Rcongas raw
saveRDS(x, paste0(out.dir, "1.rcongas_raw.rds"))
x = readRDS(paste0(out.dir, "1.rcongas_raw.rds"))

# Now plot the number of nonzero cells and the number of peaks per segment
# ggplot(x$input$segmentation, aes(x=ATAC_nonzerocells, y=ATAC_peaks)) + geom_point()
# ggplot(x$input$segmentation, aes(x=RNA_nonzerocells, y=RNA_genes)) + geom_point()

sortedATAC = x$input$segmentation %>% arrange(ATAC_nonzerocells)
sortedRNA = x$input$segmentation %>% arrange(RNA_nonzerocells)
sortedATAC
sortedRNA


segments = x$input$segmentation %>% filter(chr != 'chr13') 


x = Rcongas:::select_segments(x, segment_ids = segments$segment_id)


# Filter post mapping segments data
s_q = c(0.01, .99)
x_noout = x %>% Rcongas:::filter_outliers(lower_quantile = s_q[1], upper_quantile = s_q[2], frequency_cutoff = 0)# , frequency_cutoff = 0)#, action = 'cap')
x_back = x
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
ggsave(paste0(fig.dir, "2.1.Histogram_type.png"),
       plot_data_mod(x, to_plot = 'type', position = 'stack'),
       width = 19, height = 19)
saveRDS(x, paste0(out.dir, "2.rcongas_filtered.rds"))


#########################
ggsave(
  paste0(fig.dir, "3.1.Histogram_longest_20_segments_0outN.png"),
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
  paste0(fig.dir, "3.1.Histogram_all_segments_0out.png"),
  plot_data(
    x_noout,
    what = 'histogram',
    segments = get_input(x_noout, what = 'segmentation') %>%
      pull(segment_id)
  ),
  width = 18,
  height = 18
)

resistance = data.frame(type = unique(x_noout$input$metadata$type), resistance = c('sensitive', 'sensitive', 'resistant', 'resistant'))

x_noout$input$metadata = x$input$metadata %>% filter(cell %in% unique(x_noout$input$dataset$cell)) 


# Modificato plot_data in mdodo da plottare gli istogrammi colorati per celltype
ggsave(paste0(fig.dir, "3.1.Histogram_type_0outN.png"),
       plot_data_mod(x_noout, to_plot = 'type', position = 'stack', 
                     colors = c('red', 'yellow', 'blue', 'green', 'grey', 'orange')),
       width = 19, height = 19)


saveRDS(x_noout, paste0(out.dir, "5.rcongas_noOutliers.rds"))



