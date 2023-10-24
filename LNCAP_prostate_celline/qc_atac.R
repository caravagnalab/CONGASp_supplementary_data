
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
  # write.table(features, file = 'filtered_peak_bc_matrix/features.tsv.gz', quote = F, col.names = F, sep = '\t', row.names = F)
  
  library(Matrix)
  counts = Matrix::readMM(paste0(data_dir, '/atac/', prefix, "_matrix.mtx.gz"))
  rownames(counts) = features$peak_id
  
  
  colnames(counts) = cells$V1
  # counts <- Seurat::Read10X(data.dir = 'filtered_peak_bc_matrix', gene.column = 4)
  # metadata <- read.csv(
  #   file = "singlecell.csv",
  #   header = TRUE,
  #   row.names = 1
  # )
  #  I first ran tabix --start 1 --end 2 --zero-based
  chrom_assay <- CreateChromatinAssay(
    counts = counts,
    sep = c(":", "-"),
    genome = 'hg38',
    fragments = paste0(data_dir, '/atac/', prefix, "_fragments.tsv.gz"),
    min.cells = 10,
    min.features = 200
  )
  pbmc <- CreateSeuratObject(
    counts = chrom_assay,
    assay = "peaks"
    # meta.data = metadata
  )
  
  tmp = granges(pbmc)
  
  genome(annotations) <- "hg38"
  # add the gene information to the object
  Annotation(pbmc) <- annotations
  
  
  #QC
  # compute nucleosome signal score per cell
  pbmc <- NucleosomeSignal(object = pbmc)
  tmp = pbmc@meta.data$nucleosome_signal
  hist(tmp[tmp<5], breaks=100)
  
  # compute TSS enrichment score per cell
  pbmc <- TSSEnrichment(object = pbmc, fast = FALSE)
  
  # # add blacklist ratio and fraction of reads in peaks
  # pbmc$pct_reads_in_peaks <- pbmc$peak_region_fragments / pbmc$passed_filters * 100
  # pbmc$blacklist_ratio <- pbmc$blacklist_region_fragments / pbmc$peak_region_fragments
  
  pbmc$high.tss <- ifelse(pbmc$TSS.enrichment > 2, 'High', 'Low')
  #TSSPlot(pbmc, group.by = 'high.tss') + NoLegend()
  
  
  pbmc$nucleosome_group <- ifelse(pbmc$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
  #FragmentHistogram(object = pbmc, group.by = 'nucleosome_group')
  
  # VlnPlot(
  #   object = pbmc,
  #   features = c('nCount_peaks', 'peak_region_fragments',
  #                'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  #   pt.size = 0.1,
  #   ncol = 5
  # )
  
  pbmc_sub <- subset(
    x = pbmc,
    subset = nCount_peaks > filters_df[prefix, 'min_counts'] &
      nCount_peaks < filters_df[prefix, 'max_counts'] &
      # pct_reads_in_peaks > 15 &
      # blacklist_ratio < 0.05 &
      nucleosome_signal < 9 &
      TSS.enrichment > 2
  )
  
  pbmc = pbmc_sub
  # pbmc <- RunTFIDF(pbmc)
  # pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q0')
  # pbmc <- RunSVD(pbmc)
  
  # #DepthCor(pbmc)
  
  # pbmc <- RunUMAP(object = pbmc, reduction = 'lsi', dims = 2:30)
  # pbmc <- FindNeighbors(object = pbmc, reduction = 'lsi', dims = 2:30)
  # pbmc <- FindClusters(object = pbmc, verbose = FALSE, algorithm = 3)
  # DimPlot(object = pbmc, label = TRUE) + NoLegend()
  
  library(SummarizedExperiment)
  library(Matrix)
  se=SummarizedExperiment(assays=list(counts=pbmc@assays$peaks@counts))
  
  featureDF = as.data.frame(pbmc@assays$peaks@ranges) %>% dplyr::rename(chr = seqnames, from = start, to = end)
  
  # save_dir = paste0(prefix, '/')
  # if (!dir.exists(save_dir)){ dir.create(save_dir)}
  create_congas_tibble(sumExp=se, modality = 'ATAC', save_dir = 'atac/',
                       features = featureDF,
                       output_file = paste0(prefix, '_counts_final.tsv'))
  
  saveRDS(pbmc, paste0('atac/',prefix, '_seurat_object.rds'))
  
}


library("readxl")


segs = read_excel('CCLE_ABSOLUTE_combined_20181227.xlsx', sheet ='ABSOLUTE_combined.segtab')
head(segs)

segs = segs %>% dplyr::filter(sample == 'LNCAPCLONEFGC_PROSTATE') %>% dplyr::select(Chromosome, Start, End, Modal_Total_CN) %>%
  mutate(Chromosome = paste0('chr', Chromosome))

segs =  segs %>% dplyr::rename(chr = Chromosome, from = Start, to = End, copies = Modal_Total_CN)

write.table(segs, file = 'data/cellline_segments.tsv', quote = F, row.names = F, sep = '\t')



