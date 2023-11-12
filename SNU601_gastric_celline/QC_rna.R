machine = 'mac'


script_folder = './'
data_folder = 'data/'

source(paste0(script_folder, "/utils/congas_plot_mod.R"))
source(paste0(script_folder, "/utils/outliers_new_function.R"))
source(paste0(script_folder,'/utils/create_tibble.R'))

library(stringr)
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(patchwork)
set.seed(1234)
library(dplyr)



setwd(paste0(data_folder, "/SNU601/"))


data_dir = './'
if (!dir.exists(data_dir)){ dir.create(data_dir)}

prefixes = gsub(pattern = '_barcodes.tsv.gz', replacement = '',
                list.files(paste0(data_dir, '/rna/'), pattern = 'barcodes'))

prefix = prefixes[1] 

map_symbol_id = readr::read_delim(paste0(data_folder, '/10x_multiome/lymphoma/filtered_feature_bc_matrix/features.tsv.gz'), 
                                  delim  = '\t',
                                  col_names = c('gene_id', 'symbol', 'modality', 'chr', 'from', 'to'))  %>% dplyr::filter(modality == 'Gene Expression')


features_df = read.table(paste0(data_dir, '/rna/', prefix, "_genes.tsv.gz"))  
colnames(features_df) = c('symbol')

cells = read.table(paste0(data_dir, '/rna/', prefix, "_barcodes.tsv.gz"))

library(Matrix)

counts = Matrix::readMM(paste0(data_dir, '/rna/', prefix, "_matrix.mtx.gz"))

rownames(counts) = features_df$symbol
colnames(counts) = cells$V1

so = Seurat::CreateSeuratObject(counts = counts)

features_df = features_df %>%# dplyr::select(gene_id, gene_name) %>% 
  left_join(map_symbol_id)
duplicated_genes = features_df$symbol[duplicated(features_df$symbol)] 

counts = counts[!rownames(counts) %in% duplicated_genes,] 

features_df = features_df %>% filter(!symbol %in% duplicated_genes) %>% rename(gene=symbol)

counts = counts[features_df$gene,] 

save_dir = 'rna/'
if (!dir.exists(save_dir)){ dir.create(save_dir)}
Rcongas:::create_congas_tibble(counts=counts, modality = 'RNA', save_dir = save_dir,
                    features = features_df,
                    output_file = paste0('counts_final_nuovo.tsv'))

##### Seurat object
so <- NormalizeData(so, normalization.method = "LogNormalize", scale.factor = 10000)
all.genes <- rownames(so)
so <- ScaleData(so, features = all.genes)

so <- RunPCA(so, features = all.genes)
so <- FindNeighbors(so, dims = 1:10)
so <- FindClusters(so, resolution = 0.5)

so <- RunUMAP(so, dims = 1:10)

saveRDS(so, 'rna/seurat_object.rds')
