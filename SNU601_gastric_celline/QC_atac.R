
library(Rcongas)
library(stringi)
library(dplyr)
library(ggplot2)
library(stringr)

library(stringr)
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(patchwork)
set.seed(1234)
library(dplyr)


script_folder = './'
data_folder = 'data/'

source(paste0(script_folder, "/utils/congas_plot_mod.R"))
source(paste0(script_folder, "/utils/outliers_new_function.R"))
source(paste0(script_folder,'/utils/create_tibble.R'))

setwd(paste0(data_folder, "/SNU601"))

annotations = readRDS('annotations.rds')
data_dir = './atac/raw_data/outs/'
if (!dir.exists(data_dir)){ dir.create(data_dir)}

features = read.table(paste0(data_dir, "filtered_peak_bc_matrix/peaks.bed")) %>% 
  mutate(peak_id = paste0(V1, ':', V2, '-', V3))
cells = read.table(paste0(data_dir, "filtered_peak_bc_matrix/barcodes.tsv"))

library(Matrix)
counts = Matrix::readMM(paste0(data_dir, "filtered_peak_bc_matrix/matrix.mtx"))
rownames(counts) = features$peak_id


colnames(counts) = cells$V1
metadata <- read.csv(
  file = paste0(data_dir, "/singlecell.csv"),
  header = TRUE,
  row.names = 1
)
#  I first ran tabix --start 1 --end 2 --zero-based
chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = 'hg38',
  fragments = paste0(data_dir,  "fragments.tsv.gz"),
  min.cells = 10,
  min.features = 200
)
so <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = metadata
)

tmp = granges(so)

genome(annotations) <- "hg38"

Annotation(so) <- annotations


#QC
# compute nucleosome signal score per cell
so <- NucleosomeSignal(object = so)
tmp = so@meta.data$nucleosome_signal
hist(tmp[tmp<5], breaks=100)

# compute TSS enrichment score per cell
so <- TSSEnrichment(object = so, fast = FALSE)

# add blacklist ratio and fraction of reads in peaks
so$pct_reads_in_peaks <- so$peak_region_fragments / so$passed_filters * 100
so$blacklist_ratio <- so$blacklist_region_fragments / so$peak_region_fragments

so$high.tss <- ifelse(so$TSS.enrichment > 2, 'High', 'Low')
TSSPlot(so, group.by = 'high.tss') + NoLegend()


so$nucleosome_group <- ifelse(so$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = so, group.by = 'nucleosome_group')

VlnPlot(
  object = so_sub,
  features = c('nCount_peaks',
               'pct_reads_in_peaks',
               'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0,
  ncol = 6
)

so_sub <- subset(
  x = so,
  subset = nCount_peaks > 5000 &
    nCount_peaks < 2e5 &
    nucleosome_signal < 9 &
    TSS.enrichment > 2
)

so = so_sub


library(SummarizedExperiment)
library(Matrix)
se=SummarizedExperiment(assays=list(counts=so@assays$peaks@counts))

featureDF = as.data.frame(so@assays$peaks@ranges) %>% 
  dplyr::rename(chr = seqnames, from = start, to = end)


Rcongas:::create_congas_tibble(counts=assays(se)$counts , modality = 'ATAC', save_dir = 'atac/',
                     features = featureDF,
                     output_file = paste0('counts_final.tsv'))

saveRDS(so, paste0('atac/seurat_object.rds'))




