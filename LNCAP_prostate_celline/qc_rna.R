
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


data_dir = './'
if (!dir.exists(data_dir)){ dir.create(data_dir)}

prefixes = gsub(pattern = '_barcodes.tsv.gz', replacement = '',
                list.files(paste0(data_dir, '/rna/'), pattern = 'barcodes'))

prefix = prefixes[1] 

# I took the threshold values from the original publication
thresholds = data.frame(min_genes = c(3000, 1500, 1500, 1500), max_genes = c(7000,5000,5000,5000),
                        min_counts = c(16000, 5000, 5000, 5000), max_counts = c(50000, 25000, 25000, 25000),
                        mitoc_thresh = c(0.15,0.15,0.17,0.2))
rownames(thresholds) = prefixes

# Take feature file from 10x multiome, where each gene id is already mapped to its genomic position
map_symbol_id = readr::read_delim(paste0(data_folder, '/10x_multiome/filtered_feature_bc_matrix/features.tsv.gz'), 
                                  delim  = '\t',
                                  col_names = c('gene_id', 'symbol', 'modality', 'chr', 'from', 'to'))  %>% dplyr::filter(modality == 'Gene Expression')


for (prefix in prefixes){ 
  features_df = read.table(paste0(data_dir, '/rna/', prefix, "_features.tsv.gz"))  
  colnames(features_df) = c('gene_id', 'gene_name', 'assay')
  features_df = features_df %>% dplyr::select(gene_id, gene_name) %>% left_join(map_symbol_id, by = 'gene_id')
  
  cells = read.table(paste0(data_dir, '/rna/', prefix, "_barcodes.tsv.gz"))
  
  library(Matrix)
  
  counts = Matrix::readMM(paste0(data_dir, '/rna/', prefix, "_matrix.mtx.gz"))
  
  rownames(counts) = features_df$gene_name
  
  colnames(counts) = cells$V1
  
  so = CreateSeuratObject(counts)
  
  so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^MT-")
  so = subset(x = so, subset = nFeature_RNA > thresholds[prefix,]$min_genes &
                nFeature_RNA < thresholds[prefix,]$max_genes &
                nCount_RNA >  thresholds[prefix,]$min_counts &
                nCount_RNA < thresholds[prefix,]$max_counts &
                percent.mt < thresholds[prefix,]$mitoc_thresh * 100) 
  so@assays$RNA@meta.features  = features_df
  
  library(SummarizedExperiment)
  counts = so@assays$RNA@counts
  rownames(counts) = features_df$gene_id  
  se=SummarizedExperiment(assays=list(counts=counts))
  
  save_dir = 'rna/'
  if (!dir.exists(save_dir)){ dir.create(save_dir)}
  
  #################################### SAVE TIBBLE
  
  
  create_congas_tibble(sumExp=se, modality = 'RNA', save_dir = 'rna/',
                       features = features_df %>% dplyr::rename(gene = gene_name) ,
                       output_file = paste0(prefix, '_counts_final.tsv'))
  so = SCTransform(so, verbose = T) %>% RunPCA() %>% 
    RunUMAP(dims = 1:50)
  saveRDS(so, paste0('rna/',prefix, '_seurat_object.rds'))
  
}
