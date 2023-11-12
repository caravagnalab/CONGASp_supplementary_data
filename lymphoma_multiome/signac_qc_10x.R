
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(Signac)


script_folder = './'
data_folder = 'data/'


source('utils/create_tibble.R')

# the 10x hdf5 file contains both data types. 
h5_file = "lymph_node_lymphoma_14k_filtered_feature_bc_matrix.h5"
if (!file.exists(h5_file)) { system('curl -O https://cf.10xgenomics.com/samples/cell-arc/2.0.0/lymph_node_lymphoma_14k/lymph_node_lymphoma_14k_filtered_feature_bc_matrix.h5') }
inputdata.10x <- Read10X_h5(h5_file)

# extract RNA and ATAC data
rna_counts <- inputdata.10x$`Gene Expression`
atac_counts <- inputdata.10x$Peaks

map_symbol_id = readr::read_delim('filtered_feature_bc_matrix/features.tsv.gz', 
                                  delim  = '\t',
                                  col_names = c('gene_id', 'gene_name', 'modality', 'chr', 'from', 'to'))  %>% 
                                dplyr::filter(modality == 'Gene Expression')
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
so <- CreateSeuratObject(counts = rna_counts, meta.data = metadata)
so@assays$RNA@meta.features$gene_id = map_symbol_id$gene_id
so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^MT-")

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
so[["ATAC"]] <- chrom_assay
library(dplyr)

DefaultAssay(so) <- "ATAC"
so <- NucleosomeSignal(object = so)
tmp = so@meta.data$nucleosome_signal
hist(tmp[tmp<5], breaks=100)

# compute TSS enrichment score per cell
so <- TSSEnrichment(object = so, fast = FALSE)


so$high.tss <- ifelse(so$TSS.enrichment > 2, 'High', 'Low')
TSSPlot(so, group.by = 'high.tss') + NoLegend()

so$nucleosome_group <- ifelse(so$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = so, group.by = 'nucleosome_group')


VlnPlot(
  object = so,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5
)
DefaultAssay(so) <- "RNA"
VlnPlot(so, features = c("nCount_ATAC", "nCount_RNA","percent.mt"), ncol = 3,
        log = TRUE, pt.size = 0) + NoLegend()
so <- subset(
  x = so,
  subset = nCount_ATAC < 7e4 &
    nCount_ATAC > 5e2 &
    nCount_RNA < 25000 &
    nCount_RNA > 500 &
    percent.mt < 20 &
    TSS.enrichment > 2
)


DefaultAssay(so) <- "RNA"
so <- SCTransform(so, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

DefaultAssay(so) <- "ATAC"
so <- RunTFIDF(so)
so <- FindTopFeatures(so, min.cutoff = 'q0')
so <- RunSVD(so)
so <- RunUMAP(so, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

so <- FindMultiModalNeighbors(so, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
so <- RunUMAP(so, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
so <- FindClusters(so, graph.name = "wsnn", algorithm = 3, verbose = FALSE)

# so <- FindSubCluster(so, cluster = 6, graph.name = "wsnn", algorithm = 3)
Idents(so) <- "sub.cluster"

saveRDS(so, 'seurat_obj.rds')

so = readRDS('seurat_obj.rds')

library(SummarizedExperiment)
library(Matrix)
se=SummarizedExperiment(assays=list(counts=so@assays$ATAC@counts))

featureDF = as.data.frame(so@assays$ATAC@ranges) %>% dplyr::rename(chr = seqnames, from = start, to = end)


save_dir = 'atac/'
if (!dir.exists(save_dir)){ dir.create(save_dir)}

Rcongas:::create_congas_tibble(counts=assays(se)$counts , modality = 'ATAC', save_dir = save_dir,
                     features = featureDF,
                     output_file = 'counts_final.tsv')

se=SummarizedExperiment(assays=list(counts=so@assays$RNA@counts))

save_dir = 'rna/'
if (!dir.exists(save_dir)){ dir.create(save_dir)}


Rcongas:::create_congas_tibble(counts=assays(se)$counts , modality = 'RNA', save_dir = save_dir,
                     features = map_symbol_id %>% dplyr::rename(gene = gene_name),
                     output_file = 'counts_final.tsv')


celltypes = read.csv('lymphoma/celltypes.csv') %>% 
  dplyr::rename(cell = Barcode, type = Cell.Types )

celltypes = bind_rows(celltypes %>% mutate(cell = paste0(cell, '-ATAC')),
                      celltypes %>% mutate(cell = paste0(cell, '-RNA')))

celltypes = celltypes %>% filter(cell %in% x$input$dataset$cell)
saveRDS(celltypes, 'lymphoma/celltypes.rds')









