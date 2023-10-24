# dir.create("/mnt/storage/dati_lucreziap/r_packages_signac/")
# .libPaths("/mnt/storage/dati_lucreziap/r_packages_signac/")

# install.packages("BiocManager")

BiocManager::install('GenomeInfoDb')
setwd("/mnt/storage/dati_lucreziap/CONGASp_data/copy_scAT/adult_4218/adult4218/outs")

# setRepositories(ind=1:3) # needed to automatically install Bioconductor dependencies
# install.packages("SeuratObject")
# BiocManager::install(c('BSgenome.Hsapiens.UCSC.hg38', 'EnsDb.Hsapiens.v86'))
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(patchwork)
set.seed(1234)
library(dplyr)

features = read.table("filtered_peak_bc_matrix/peaks.bed") %>% 
  mutate(peak_id = paste0(V1, ':', V2, '-', V3))

# write.table(features, file = 'filtered_peak_bc_matrix/features.tsv.gz', quote = F, col.names = F, sep = '\t', row.names = F)

library(Matrix)
counts = Matrix::readMM('filtered_peak_bc_matrix/matrix.mtx')
rownames(counts) = features$peak_id

cells = read.table('filtered_peak_bc_matrix/barcodes.tsv')
colnames(counts) = cells$V1
# counts <- Seurat::Read10X(data.dir = 'filtered_peak_bc_matrix', gene.column = 4)
metadata <- read.csv(
  file = "singlecell.csv",
  header = TRUE,
  row.names = 1
)

chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = 'hg38',
  fragments = 'fragments.tsv.gz',
  min.cells = 10,
  min.features = 200
)
pbmc <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = metadata
)

# To see the genomic ranges associated with each feature in the object:
tmp = granges(pbmc)
# Add gene annotations to the pbmc object for the human genome:
# annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
# # change to UCSC style since the data was mapped to hg19
# seqlevelsStyle(annotations) <- 'UCSC'
annotations = readRDS('/mnt/storage/dati_lucreziap/CONGASp_data/10x_multiome/lymphoma/annotations.rds')
genome(annotations) <- "hg38"

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

# add blacklist ratio and fraction of reads in peaks
pbmc$pct_reads_in_peaks <- pbmc$peak_region_fragments / pbmc$passed_filters * 100
pbmc$blacklist_ratio <- pbmc$blacklist_region_fragments / pbmc$peak_region_fragments

pbmc$high.tss <- ifelse(pbmc$TSS.enrichment > 2, 'High', 'Low')
TSSPlot(pbmc, group.by = 'high.tss') + NoLegend()


pbmc$nucleosome_group <- ifelse(pbmc$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = pbmc, group.by = 'nucleosome_group')

VlnPlot(
  object = pbmc,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5
)

pbmc_sub <- subset(
  x = pbmc,
  subset = peak_region_fragments > 3000 &
    peak_region_fragments < 50000 &
    pct_reads_in_peaks > 15 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2
)
pbmc_sub

pbmc = pbmc_sub
pbmc <- RunTFIDF(pbmc)
pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q0')
pbmc <- RunSVD(pbmc)

DepthCor(pbmc)

pbmc <- RunUMAP(object = pbmc, reduction = 'lsi', dims = 2:30)
pbmc <- FindNeighbors(object = pbmc, reduction = 'lsi', dims = 2:30)
pbmc <- FindClusters(object = pbmc, verbose = FALSE, algorithm = 3)
DimPlot(object = pbmc, label = TRUE) + NoLegend()

saveRDS(pbmc, '/mnt/storage/dati_lucreziap/CONGASp_data/copy_scAT/adult_4218/seurat_object.rds')





