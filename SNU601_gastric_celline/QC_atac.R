# dir.create("/mnt/storage/dati_lucreziap/r_packages_signac/")
# .libPaths("/mnt/storage/dati_lucreziap/r_packages_signac/")

# install.packages("BiocManager")

# BiocManager::install('GenomeInfoDb')
# setwd("/mnt/storage/dati_lucreziap/CONGASp_data/copy_scAT/adult_4218/adult4218/outs")

# setRepositories(ind=1:3) # needed to automatically install Bioconductor dependencies
# install.packages("SeuratObject")
# BiocManager::install(c('BSgenome.Hsapiens.UCSC.hg38', 'EnsDb.Hsapiens.v86'))
machine = 'guascone'
if (machine == 'guascone')
  dyn.load('/mnt/storage/dati_lucreziap/conda_envs/cn_atac/lib/libiconv.so.2')
library(Rcongas)
library(stringi)
library(dplyr)
library(ggplot2)
library(stringr)

if (machine == 'guascone'){ 
  script_folder = '/mnt/storage/dati_lucreziap/congasp_analysis/'
  data_folder = '/mnt/storage/dati_lucreziap/CONGASp_data/'
} else if (machine == 'giacinto'){ 
  script_folder = '/dati-raid/BimiB/lpatruno/ATAC/congasp_analysis'
  data_folder = '/dati-raid/BimiB/lpatruno/ATAC/CONGASp_data/'
} else if (machine == 'mac'){
  script_folder = '/Users/lucrezia/Library/CloudStorage/OneDrive-UniversityCollegeLondon/dottorato/congasp_analysis/'
  data_folder = '/Users/lucrezia/Dropbox/CONGASp_data/'
}

library(stringr)
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(patchwork)
set.seed(1234)
library(dplyr)

if (machine == 'guascone'){ 
  script_folder = '/mnt/storage/dati_lucreziap/congasp_analysis/'
  data_folder = '/mnt/storage/dati_lucreziap/CONGASp_data/'
} else if (machine == 'giacinto'){ 
  script_folder = '/dati-raid/BimiB/lpatruno/ATAC/congasp_analysis'
  data_folder = '/dati-raid/BimiB/lpatruno/ATAC/CONGASp_data/'
}

source(paste0(script_folder, "/congas_plot_mod.R"))
source(paste0(script_folder, "/outliers_new_function.R"))
source(paste0(script_folder,'/analyse_data/create_tibble.R'))

setwd(paste0(data_folder, "/Reviews/SNU650"))

# Add gene annotations to the pbmc object for the human genome:
#annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
# # change to UCSC style since the data was mapped to hg19
#seqlevelsStyle(annotations) <- 'UCSC'
annotations = readRDS(paste0(data_folder,'10x_multiome/lymphoma/annotations.rds'))
data_dir = './atac/raw_data/outs/'
if (!dir.exists(data_dir)){ dir.create(data_dir)}

features = read.table(paste0(data_dir, "filtered_peak_bc_matrix/peaks.bed")) %>% 
  mutate(peak_id = paste0(V1, ':', V2, '-', V3))
cells = read.table(paste0(data_dir, "filtered_peak_bc_matrix/barcodes.tsv"))
# write.table(features, file = 'filtered_peak_bc_matrix/features.tsv.gz', quote = F, col.names = F, sep = '\t', row.names = F)

library(Matrix)
counts = Matrix::readMM(paste0(data_dir, "filtered_peak_bc_matrix/matrix.mtx"))
rownames(counts) = features$peak_id


colnames(counts) = cells$V1
# counts <- Seurat::Read10X(data.dir = 'filtered_peak_bc_matrix', gene.column = 4)
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
pbmc <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = metadata
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

# add blacklist ratio and fraction of reads in peaks
pbmc$pct_reads_in_peaks <- pbmc$peak_region_fragments / pbmc$passed_filters * 100
pbmc$blacklist_ratio <- pbmc$blacklist_region_fragments / pbmc$peak_region_fragments

pbmc$high.tss <- ifelse(pbmc$TSS.enrichment > 2, 'High', 'Low')
TSSPlot(pbmc, group.by = 'high.tss') + NoLegend()


pbmc$nucleosome_group <- ifelse(pbmc$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = pbmc, group.by = 'nucleosome_group')

VlnPlot(
  object = pbmc_sub,
  features = c('nCount_peaks',
               'pct_reads_in_peaks',
               'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0,
  ncol = 6
)

pbmc_sub <- subset(
  x = pbmc,
  subset = nCount_peaks > 5000 &
    nCount_peaks < 2e5 &
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

featureDF = as.data.frame(pbmc@assays$peaks@ranges) %>% 
  dplyr::rename(chr = seqnames, from = start, to = end)

# save_dir = paste0(prefix, '/')
# if (!dir.exists(save_dir)){ dir.create(save_dir)}
Rcongas:::create_congas_tibble(counts=assays(se)$counts , modality = 'ATAC', save_dir = 'atac/',
                     features = featureDF,
                     output_file = paste0('counts_final.tsv'))

saveRDS(pbmc, paste0('atac/seurat_object.rds'))




