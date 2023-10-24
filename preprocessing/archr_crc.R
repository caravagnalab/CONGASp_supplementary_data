dyn.load('/mnt/storage/dati_lucreziap/conda_envs/cn_atac/lib/libgeos_c.so.1')
dyn.load('/mnt/storage/dati_lucreziap/conda_envs/cn_atac/lib/rgeos.so.1') 
# Sys.unsetenv('GITHUB_PAT')
# devtools::install_github("GreenleafLab/ArchR", ref="master", repos = BiocManager::repositories())
# options(timeout = 1000)
# BiocManager::install('BSgenome.Hsapiens.UCSC.hg38')
setwd('/mnt/storage/dati_lucreziap/CONGASp_data/crc')
source('/mnt/storage/dati_lucreziap/congasp_script/analyse_data/create_tibble.R')
library(BSgenome.Hsapiens.UCSC.hg38)
library(ArchR)

addArchRGenome("hg38")

patients = c('A001', 'A015')

conditions = c('normal', 'polyp')#'tumor' 
patient = 'A001'
tumorNormal = 'tumor'

output_dir = 'arrow/QualityControl'
if (!dir.exists(output_dir)) { dir.create(output_dir)}

# Thresholds for mifragments:
#  A001: normal = 10^3.5, tumor = 10^3, polyp = 10^3.5

for (patient in patients){ 
  for (tumorNormal in conditions){
    
    save_dir = paste0(patient, '/', tumorNormal, '/atac')
    # if (!file.exists(paste0(save_dir, '/counts_final.tsv'))){ 
    ArrowFiles <- createArrowFiles(
      inputFiles = paste0(patient, '/', tumorNormal, '/atac/fragments.tsv.gz'),
      sampleNames = paste(patient, tumorNormal, sep = '_'),
      minTSS = 6, #Dont set this too high because you can always increase later
      filterFrags = 5000, 
      QCDir = output_dir,
      force = F,
      addTileMat = TRUE,
      addGeneScoreMat = FALSE,
      excludeChr = c('chrM'),
      TileMatParams = list('binarize' = FALSE)
    )
    
    # prova = readRDS('/mnt/storage/dati_lucreziap/CONGASp_data/crc/arrow/QualityControl/A001_tumor/A001_tumor-Pre-Filter-Metadata.rds')
    # 
    # prova = as.data.frame(prova)
    
    # doubScores <- addDoubletScores(
    #   input = ArrowFiles,
    #   k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
    #   knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search.
    #   LSIMethod = 1,
    #   outDir = output_dir
    # )
    proj <- ArchRProject(
      ArrowFiles = ArrowFiles, 
      outputDirectory = output_dir,
      copyArrows = TRUE #This is recommened so that you maintain an unaltered copy for later usage.
    )
    
    prova = h5read(ArrowFiles, '/TileMatrix/')
    
    proj <- addIterativeLSI(ArchRProj = proj, useMatrix = "TileMatrix", name = "IterativeLSI")
    proj <- addClusters(input = proj, reducedDims = "IterativeLSI")
    proj <- addUMAP(ArchRProj = proj, reducedDims = "IterativeLSI")
    p2 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
    # p2
    
    plotPDF(p2, name = paste0(patient, '_', tumorNormal, "_Plot-UMAP-Clusters.pdf"),
            ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)
    
    
    # counts = as.matrix(assay(sim))
    
    # featureDF <- ArchR:::.getFeatureDF(ArrowFiles, subGroup = 'TileMatrix')
    # seqnames <- ArchR:::.availableSeqnames(ArrowFiles, subGroup = 'TileMatrix')
    
    
    featureDF = tibble::tibble(h5read(ArrowFiles, '/TileMatrix/Info/FeatureDF')) %>% mutate(end = start + 500) %>%
      mutate(start = start + 1) %>% mutate(tileID = paste0(seqnames, ':', start, '-', end))
    featureDF = featureDF %>% dplyr::rename(chr = seqnames, from = start, to = end)
    sim = getMatrixFromProject(proj,#  useSeqnames = 'chr11', 
                               useMatrix = 'TileMatrix')
    
    rownames(sim) = featureDF$tileID
    
    create_congas_tibble(sumExp=sim, modality = 'ATAC', save_dir = save_dir,
                         features = featureDF,
                         output_file = 'counts_final.tsv')
    # } else { 
    #   print(paste0(patient, ' ', tumorNormal, ' already found'))
    #   }
  } 
}

patient = 'A001'
tumorNormal = 'normal'

save_dir = paste0(patient, '/', tumorNormal, '/atac2')
if (!dir.exists(save_dir)){ dir.create(save_dir)}
# if (!file.exists(paste0(save_dir, '/counts_final.tsv'))){ 
ArrowFiles <- createArrowFiles(
  inputFiles = c(paste0(patient, '/', tumorNormal, '/atac/fragments.tsv.gz'),
                 paste0(patient, '/', tumorNormal, '/atac/fragments2.tsv.gz')),
  sampleNames = c(paste(patient, tumorNormal, sep = '_'),
                  paste(patient, tumorNormal, '2', sep = '_')),
  minTSS = 4, #Dont set this too high because you can always increase later
  filterFrags = 2500, 
  QCDir = output_dir,
  force = T,
  addTileMat = TRUE,
  addGeneScoreMat = FALSE,
  excludeChr = c('chrM'),
  TileMatParams = list('binarize' = FALSE)
)


proj <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = output_dir,
  copyArrows = TRUE #This is recommened so that you maintain an unaltered copy for later usage.
)

prova = h5read(ArrowFiles, '/TileMatrix/')

proj <- addIterativeLSI(ArchRProj = proj, useMatrix = "TileMatrix", name = "IterativeLSI")
proj <- addClusters(input = proj, reducedDims = "IterativeLSI")
proj <- addUMAP(ArchRProj = proj, reducedDims = "IterativeLSI")
p2 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
# p2

plotPDF(p2, name = paste0(patient, '_', tumorNormal, "_Plot-UMAP-Clusters.pdf"),
        ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)
p3 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Sample", embedding = "UMAP")

plotPDF(p3, name = paste0(patient, '_', tumorNormal, "_Plot-UMAP-Sample.pdf"),
        ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

featureDF = tibble::tibble(h5read(ArrowFiles[1] , '/TileMatrix/Info/FeatureDF')) %>% mutate(end = start + 500) %>%
  mutate(start = start + 1) %>% mutate(tileID = paste0(seqnames, ':', start, '-', end))

# featureDF2 = tibble::tibble(h5read(ArrowFiles[2] , '/TileMatrix/Info/FeatureDF')) %>% mutate(end = start + 500) %>%
#   mutate(start = start + 1) %>% mutate(tileID = paste0(seqnames, ':', start, '-', end))

featureDF = featureDF %>% dplyr::rename(chr = seqnames, from = start, to = end)
sim = getMatrixFromProject(proj,#  useSeqnames = 'chr11', 
                           useMatrix = 'TileMatrix')

rownames(sim) = featureDF$tileID

create_congas_tibble(sumExp=sim, modality = 'ATAC', save_dir = save_dir,
                     features = featureDF,
                     output_file = 'counts_final.tsv')

