# BiocManager::install('Signac')
# install.packages('hdf5r', configure.args="--with-hdf5=/mnt/storage/dati_lucreziap/conda_envs/cn_atac/bin/h5cc")
# BiocManager::install("EnsDb.Hsapiens.v86")
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(Signac)
setwd('/mnt/storage/dati_lucreziap/CONGASp_data/10x_multiome/lymphoma')
dyn.load('/mnt/storage/dati_lucreziap/conda_envs/cn_atac/lib/libiconv.so.2')
dyn.load('/mnt/storage/dati_lucreziap/conda_envs/cn_atac/lib/libhdf5_hl.so.100')
# Sys.setenv('PKG_CONFIG_PATH' = '/mnt/storage/dati_lucreziap/conda_envs/cn_atac/lib/pkgconfig/')
Sys.setenv('LD_LIBRARY_PATH' = paste0(Sys.getenv('LD_LIBRARY_PATH'),
                                      ':/mnt/storage/dati_lucreziap/conda_envs/cn_atac/lib'))
# the 10x hdf5 file contains both data types. 
h5_file = "lymph_node_lymphoma_14k_filtered_feature_bc_matrix.h5"
if (!file.exists(h5_file)) { system('curl -O https://cf.10xgenomics.com/samples/cell-arc/2.0.0/lymph_node_lymphoma_14k/lymph_node_lymphoma_14k_filtered_feature_bc_matrix.h5') }
inputdata.10x <- Read10X_h5(h5_file)

# extract RNA and ATAC data
rna_counts <- inputdata.10x$`Gene Expression`
atac_counts <- inputdata.10x$Peaks

map_symbol_id = readr::read_delim('filtered_feature_bc_matrix/features.tsv.gz', delim  = '\t',
                                  col_names = c('gene_id', 'gene_name', 'modality', 'chr', 'from', 'to'))  %>% dplyr::filter(modality == 'Gene Expression')
dim(map_symbol_id)
dim(rna_counts)
rownames(rna_counts)[map_symbol_id$gene_name != rownames(rna_counts)]

# Read metadata (I added this, to compute TSS enrichment and nucleosome signal)
metadata <- read.csv(
  file = "lymph_node_lymphoma_14k_per_barcode_metrics.csv",
  header = TRUE,
  row.names = 1
) %>% dplyr::filter(gex_barcode %in% colnames(rna_counts))
# Create Seurat object
pbmc <- CreateSeuratObject(counts = rna_counts, meta.data = metadata)
pbmc@assays$RNA@meta.features$gene_id = map_symbol_id$gene_id
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Now add in the ATAC-seq data
# we'll only use peaks in standard chromosomes
# grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
# grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
# atac_counts <- atac_counts[as.vector(grange.use), ]
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
ensembldb::seqlevelsStyle(annotations) <- 'UCSC'
annotations = readRDS('annotations.rds')
genome(annotations) <- "hg38"

frag.file <- "lymph_node_lymphoma_14k_atac_fragments.tsv.gz"
if (!file.exists(frag.file)) { 
  system('curl -O https://cf.10xgenomics.com/samples/cell-arc/2.0.0/lymph_node_lymphoma_14k/lymph_node_lymphoma_14k_atac_fragments.tsv.gz') 
  system('curl -O https://cf.10xgenomics.com/samples/cell-arc/2.0.0/lymph_node_lymphoma_14k/lymph_node_lymphoma_14k_atac_fragments.tsv.gz.tbi')
}
chrom_assay <- CreateChromatinAssay(
  counts = atac_counts,
  sep = c(":", "-"),
  genome = 'hg38',
  fragments = frag.file,
  min.cells = 10,
  annotation = annotations
)
pbmc[["ATAC"]] <- chrom_assay
# We perform basic QC based on the num


# We next perform pre-processing and dimensional reduction on both assays independently, using standard approaches for RNA and ATAC-seq data.
library(dplyr)
# RNA analysis

# ATAC analysis
# We exclude the first dimension as this is typically correlated with sequencing depth
DefaultAssay(pbmc) <- "ATAC"
pbmc <- NucleosomeSignal(object = pbmc)
tmp = pbmc@meta.data$nucleosome_signal
hist(tmp[tmp<5], breaks=100)

# compute TSS enrichment score per cell
pbmc <- TSSEnrichment(object = pbmc, fast = FALSE)

# add blacklist ratio and fraction of reads in peaks
# pbmc$pct_reads_in_peaks <- pbmc$peak_region_fragments / pbmc$passed_filters * 100
# pbmc$blacklist_ratio <- pbmc$blacklist_region_fragments / pbmc$peak_region_fragments

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
DefaultAssay(pbmc) <- "RNA"
VlnPlot(pbmc, features = c("nCount_ATAC", "nCount_RNA","percent.mt"), ncol = 3,
        log = TRUE, pt.size = 0) + NoLegend()
pbmc <- subset(
  x = pbmc,
  subset = nCount_ATAC < 7e4 &
    nCount_ATAC > 5e2 &
    nCount_RNA < 25000 &
    nCount_RNA > 500 &
    percent.mt < 20 &
    TSS.enrichment > 2
)
####### Aggiunta per salvare oggetto tutorial
rna_cells = readRDS('tutorial/congas_data/rna_cells.rds') %>% 
  gsub(pattern = '-RNA', replacement = '', x = .)
atac_cells = readRDS('tutorial/congas_data/atac_cells.rds') %>% 
  gsub(pattern = '-ATAC', replacement = '', x = .)
pbmc_sub = pbmc[,rna_cells] 

rna_counts = pbmc_sub@assays$RNA@counts  
atac_counts = pbmc_sub@assays$ATAC@counts  

saveRDS(rna_counts, 'tutorial/rna_counts.rds')
saveRDS(atac_counts, 'tutorial/atac_counts.rds')


#######


DefaultAssay(pbmc) <- "RNA"
pbmc <- SCTransform(pbmc, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

DefaultAssay(pbmc) <- "ATAC"
pbmc <- RunTFIDF(pbmc)
pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q0')
pbmc <- RunSVD(pbmc)
pbmc <- RunUMAP(pbmc, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

pbmc <- FindMultiModalNeighbors(pbmc, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
pbmc <- RunUMAP(pbmc, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
pbmc <- FindClusters(pbmc, graph.name = "wsnn", algorithm = 3, verbose = FALSE)

# pbmc <- FindSubCluster(pbmc, cluster = 6, graph.name = "wsnn", algorithm = 3)
Idents(pbmc) <- "sub.cluster"
library(ggplot2)

p1 <- DimPlot(pbmc, reduction = "umap.rna", group.by = "seurat_clusters", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(pbmc, reduction = "umap.atac", group.by = "seurat_clusters", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(pbmc, reduction = "wnn.umap", group.by = "seurat_clusters", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
ggsave(plot = p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5)), filename = 'umap_clusters.png',
       width = 10, height = 5)

ggsave(plot = p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5)), filename = 'umap_clusters.pdf',
       width = 10, height = 5)

saveRDS(pbmc, 'seurat_obj.rds')

so = readRDS('seurat_obj.rds')

library(SummarizedExperiment)
library(Matrix)
se=SummarizedExperiment(assays=list(counts=pbmc@assays$ATAC@counts))

featureDF = as.data.frame(pbmc@assays$ATAC@ranges) %>% dplyr::rename(chr = seqnames, from = start, to = end)


save_dir = 'atac/'
if (!dir.exists(save_dir)){ dir.create(save_dir)}

#################################### SAVE TIBBLE


source('/mnt/storage/dati_lucreziap/congasp_analysis/analyse_data/create_tibble.R')
Rcongas:::create_congas_tibble(counts=assays(se)$counts , modality = 'ATAC', save_dir = save_dir,
                     features = featureDF,
                     output_file = 'counts_final.tsv')

se=SummarizedExperiment(assays=list(counts=pbmc@assays$RNA@counts))

save_dir = 'rna/'
if (!dir.exists(save_dir)){ dir.create(save_dir)}

#################################### SAVE TIBBLE


Rcongas:::create_congas_tibble(counts=assays(se)$counts , modality = 'RNA', save_dir = save_dir,
                     features = map_symbol_id %>% dplyr::rename(gene = gene_name),
                     output_file = 'counts_final.tsv')


celltypes = read.csv('lymphoma/tutorial/celltypes.csv') %>% 
  dplyr::rename(cell = Barcode, type = Cell.Types )

celltypes = bind_rows(celltypes %>% mutate(cell = paste0(cell, '-ATAC')),
                      celltypes %>% mutate(cell = paste0(cell, '-RNA')))

celltypes = celltypes %>% filter(cell %in% x$input$dataset$cell)
saveRDS(celltypes, 'lymphoma/tutorial/celltypes.rds')









