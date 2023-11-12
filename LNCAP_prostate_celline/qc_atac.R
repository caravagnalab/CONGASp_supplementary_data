
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

setwd(paste0(data_folder, "/celline_prostate"))

# Add gene annotations to the pbmc object for the human genome:

annotations = readRDS('annotations.rds')
data_dir = './'
if (!dir.exists(data_dir)){ dir.create(data_dir)}

prefixes = gsub(pattern = '_barcodes.tsv.gz', replacement = '',
                list.files(paste0(data_dir, '/atac/'), pattern = 'barcodes'))

filters_df = data.frame(pct_reads_peaks = c(30, 30, 40, 30),
                        min_counts = c(2000, 1000, 1000, 1000),
                        max_counts = c(20000, 20000, 20000, 25000)
                         )
rownames(filters_df) = prefixes
prefix = prefixes[4]
for (prefix in prefixes) { 
  features = read.table(paste0(data_dir, '/atac/', prefix, "_peaks.bed.gz")) %>% 
    mutate(peak_id = paste0(V1, ':', V2, '-', V3))
  cells = read.table(paste0(data_dir, '/atac/', prefix, "_barcodes.tsv.gz"))
  
  
  library(Matrix)
  counts = Matrix::readMM(paste0(data_dir, '/atac/', prefix, "_matrix.mtx.gz"))
  rownames(counts) = features$peak_id
  
  
  colnames(counts) = cells$V1
  
  chrom_assay <- CreateChromatinAssay(
    counts = counts,
    sep = c(":", "-"),
    genome = 'hg38',
    fragments = paste0(data_dir, '/atac/', prefix, "_fragments.tsv.gz"),
    min.cells = 10,
    min.features = 200
  )
  so <- CreateSeuratObject(
    counts = chrom_assay,
    assay = "peaks"
  )
  
  tmp = granges(so)
  
  genome(annotations) <- "hg38"
  # add the gene information to the object
  Annotation(so) <- annotations
  
  
  #QC
  # compute nucleosome signal score per cell
  so <- NucleosomeSignal(object = so)
  tmp = so@meta.data$nucleosome_signal
  hist(tmp[tmp<5], breaks=100)
  
  # compute TSS enrichment score per cell
  so <- TSSEnrichment(object = so, fast = FALSE)
  
  
  so$high.tss <- ifelse(so$TSS.enrichment > 2, 'High', 'Low')
  
  
  so$nucleosome_group <- ifelse(so$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
  
  
  so_sub <- subset(
    x = so,
    subset = nCount_peaks > filters_df[prefix, 'min_counts'] &
      nCount_peaks < filters_df[prefix, 'max_counts'] &
      nucleosome_signal < 9 &
      TSS.enrichment > 2
  )
  
  so = so_sub
  
  
  library(SummarizedExperiment)
  library(Matrix)
  se=SummarizedExperiment(assays=list(counts=so@assays$peaks@counts))
  
  featureDF = as.data.frame(so@assays$peaks@ranges) %>% dplyr::rename(chr = seqnames, from = start, to = end)
  
  create_congas_tibble(sumExp=se, modality = 'ATAC', save_dir = 'atac/',
                       features = featureDF,
                       output_file = paste0(prefix, '_counts_final.tsv'))
  
  saveRDS(so, paste0('atac/',prefix, '_seurat_object.rds'))
  
}


library("readxl")


segs = read_excel('CCLE_ABSOLUTE_combined_20181227.xlsx', sheet ='ABSOLUTE_combined.segtab')
head(segs)

segs = segs %>% dplyr::filter(sample == 'LNCAPCLONEFGC_PROSTATE') %>% dplyr::select(Chromosome, Start, End, Modal_Total_CN) %>%
  mutate(Chromosome = paste0('chr', Chromosome))

segs =  segs %>% dplyr::rename(chr = Chromosome, from = Start, to = End, copies = Modal_Total_CN)

write.table(segs, file = 'data/cellline_segments.tsv', quote = F, row.names = F, sep = '\t')



