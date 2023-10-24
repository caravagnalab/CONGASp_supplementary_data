machine = 'mac'

if (machine == 'guascone') {
  dyn.load('/mnt/storage/dati_lucreziap/conda_envs/cn_atac/lib/libiconv.so.2')
  library(reticulate)
  reticulate::use_condaenv('/mnt/storage/dati_lucreziap/conda_envs/r-reticulate')
  script_folder = '/mnt/storage/dati_lucreziap/congasp_analysis/'
  data_folder = '/mnt/storage/dati_lucreziap/CONGASp_data/'
} else if (machine == 'giacinto') {
  library(reticulate)
  reticulate::use_condaenv('/dati-raid/BimiB/share/conda_envs/congasp_reticulate/')
  script_folder = '/dati-raid/BimiB/lpatruno/ATAC/congasp_analysis'
  data_folder = '/dati-raid/BimiB/lpatruno/ATAC/CONGASp_data/'
} else if (machine == 'docker') {
  library(reticulate)
  reticulate::use_condaenv('/home/root/anaconda3')
  script_folder = '/app/congasp_analysis/'
  data_folder = '/app/'
} else if (machine == 'mac'){
  script_folder = '/Users/lucrezia/Library/CloudStorage/OneDrive-UniversityCollegeLondon/dottorato/congasp_analysis/'
  data_folder = '/Users/lucrezia/Dropbox/CONGASp_data/'
}
library(Rcongas)
library(stringi)
library(dplyr)
library(ggplot2)
library(stringr)
library(doParallel)
library(tidyverse)
source(paste0(script_folder, "/congas_plot_mod.R"))
source(paste0(script_folder, "/plot_cna_clusters.R"))


setwd(paste0(data_folder, "/Reviews/SNU601/"))

lik_rna = 'NB'
lik_atac = 'NB'
prefix = paste0('bimodal_rna_', lik_rna, '_atac_', lik_atac, "/")
library(tidyr)

# obj = readRDS(paste0(prefix, "congas_data/4.rcongas_noOutliers.rds"))
# obj$input$normalisation = rbind(auto_normalisation_factor(obj$input$dataset %>% filter(modality == 'ATAC')) %>% mutate(modality = 'ATAC'),
# 					auto_normalisation_factor(obj$input$dataset %>% filter(modality == 'RNA')) %>% mutate(modality = 'RNA'))
# celltypes = readRDS('bimodal_rna_NB_atac_NB/congas_data/0.celltypes.rds')
#  obj$input$metadata = celltypes                                                     
# p = plot_data_histogram_mod(obj, to_plot = 'type', position = 'stack')
#   ggsave(paste0(out.dir, '/prova_renorm.pdf'), p, height=15, width=15)  


x_filt = readRDS(paste0(prefix, "congas_data/5.rcongas_noOutliers.rds"))

cellularity = NULL#0.7


k = c(1:8)

binom_limits = c(40,1000)
model = "BIC"
steps_try = c(15000)

lr = 0.01
temperature = 20#10

steps = steps_try[1]
lambda = 0.5#c(0.2,0.5,0.7) #seq(0.05, 0.95, 0.1)

x_filt = Rcongas:::select_segments(x_filt, 
                                  segment_ids = setdiff(x_filt$input$segmentation$segment_id, 
                                                     'chr13:80640000:95100000')) 

out_folder_prefix = paste0('NEWASScongas_run_rna_', paste(k, collapse = '-'),'_lambda_', lambda, '_sameMix/')
out_folder = paste0(prefix, "/congas_figures/", out_folder_prefix)
out_folder_data = paste0(prefix, "/congas_data/", out_folder_prefix)

if (!dir.exists(out_folder)) {dir.create(out_folder, recursive = T)}
if (!dir.exists(out_folder_data)) {dir.create(out_folder_data, recursive = T)}

hyperparams_filt <- auto_config_run(x_filt, k, 
                                    prior_cn=c(0.2, 0.6, 0.1, 0.05, 0.05),
                                    purity = cellularity, init_importance = 0.6)

hyperparams_filt$b = 1
hyperparams_filt$a = 0.1
hyperparams_filt$binom_prior_limits = binom_limits

fit_filt <- Rcongas:::fit_congas(x_filt,
                                 K = k, 
                                 lambdas = lambda, 
                                 learning_rate = lr, 
                                 steps = steps,
                                 model_parameters = hyperparams_filt, 
                                 model_selection = model,
                                 latent_variables = "G",
                                 compile = F,
                                 temperature = temperature, 
                                 same_mixing = T, 
                                 threshold = 0.001)


saveRDS(fit_filt, file=paste0(out_folder_data, "fit.rds"))
fit_filt = readRDS(paste0(out_folder_data, "fit.rds"))
#pdf(paste0(out_folder, "density.pdf"))
pdf(paste0(out_folder, "summary.pdf"))
p = Rcongas::plot_fit(fit_filt, highlights = T, what='density')
print(p)
p = Rcongas::plot_fit(fit_filt, what='posterior_CNA')
print(p)
p = Rcongas::plot_fit(fit_filt, what='heatmap')
print(p)
p = Rcongas::plot_fit(fit_filt, what='CNA')
print(p)
p = Rcongas::plot_fit(fit_filt, what='scores')
print(p)
dev.off()

pdf(paste0(out_folder, "summary_densities.pdf"))
p = plot_fit(fit_filt, highlights = F, what = 'density')
print(p)
dev.off()


offsets = tibble(clone = c(1:6), offset = seq(0.1, 0.6, by = 0.1))

# Without chr13. clone 3 and clone 6 are the same. So, I remove clone 6 from the ground truth

genotypes = readRDS('atac/cloneGenotypes.rds') %>% 
  mutate(segment_id = paste0(chr, ':', from, ':', to)) %>% 
  select(copies, clone, segment_id) %>% left_join(offsets) %>%
  mutate(copies_plot = copies + offset) %>% 
  filter(clone != 6)

gen_plot = genotypes %>% 
  mutate(clone = as.factor(clone)) %>% separate(segment_id, into = c('chr', 'from', 'to'))

p = CNAqc:::blank_genome() +
  geom_rect(data = gen_plot ,
            mapping = aes(xmin = from, xmax=to, ymin = copies, ymax = copies+0.1 ,
                          color = clone)) 
p
segmentation = fit_filt$input$segmentation
#Associate to every inferred cluster the closest one in the ground truth based on inferred copy number values: 
clone_cluster = fit_filt$best_fit$CNA %>% 
  left_join(genotypes) %>% 
  group_by(cluster, clone) %>% 
  summarise(difference = sum(abs(value-copies)) / nrow(segmentation)) %>% 
  filter(difference == min(difference))

inferred_gt_cn = fit_filt$best_fit$CNA %>% 
  left_join(clone_cluster) %>% 
  left_join(genotypes) %>%  filter(! cluster %in% c('C3', 'C6')) %>%
  separate(segment_id, into = c('chr','from','to'), sep = ':') %>%
  mutate(from = as.integer(from), to = as.integer(to)) %>% 
  dplyr::rename(CONGASp = value, scDNA = copies_plot) %>%
  pivot_longer(cols = c(CONGASp, scDNA), names_to = 'Source')

reference_genome = CNAqc:::get_reference('GRCh38')

vfrom = reference_genome$from
names(vfrom) = reference_genome$chr

inferred_gt_cn = inferred_gt_cn %>% 
  mutate(from = from + vfrom[chr],
         to = to + vfrom[chr])

p = CNAqc:::blank_genome() +
  geom_rect(data = inferred_gt_cn %>% mutate(clone = as.factor(clone)),
            mapping = aes(xmin = from, xmax=to, ymin = value, ymax = value+0.1, color = clone, fill = Source)) +
  facet_wrap(~cluster)

ggsave(filename = paste0(out_folder, 'inferred_vs_gt.pdf'), plot = p)

gt_genotypes = genotypes %>% 
  separate(segment_id, into = c('chr','from','to'), sep = ':', remove = F) %>%
  mutate(from = as.integer(from), to = as.integer(to), clone = as.character(clone))  %>%
  bind_rows(fit_filt$best_fit$CNA %>% 
              separate(segment_id, into = c('chr','from','to'), sep = ':', remove = F) %>%
              mutate(from = as.integer(from), to = as.integer(to)) %>%
              rename(clone = cluster, copies = value))

reference_genome = CNAqc:::get_reference('GRCh38')

vfrom = reference_genome$from
names(vfrom) = reference_genome$chr

gt_genotypes = gt_genotypes %>% 
  mutate(from = from + vfrom[chr],
         to = to + vfrom[chr],
         copies = as.factor(copies))

ggplot(gt_genotypes, aes(x = segment_id, y = clone, fill = copies, label = copies)) + 
  geom_tile() + geom_text() + scale_fill_brewer(palette="Dark2")




bms_idx <- order(fit_filt$model_selection %>%  pull('BIC'))
fit_filt$model_selection <-  fit_filt$model_selection[bms_idx,]

myCluster <- makeCluster(10, 
                         type = "FORK")
registerDoParallel(myCluster)

# foreach (r = seq(1, nrow(fit_filt$model_selection))) %dopar% {
for (r in seq(1, nrow(fit_filt$model_selection))) {
  #i = fit_filt$model_selection[r,]$run_index
  i = r
  print(i)
  k = fit_filt$model_selection$hyperparameter_K[r]
  lambda = fit_filt$model_selection$lambda[r]
  best_fit  <-  Rcongas:::format_best_model(fit_filt, fit_filt$runs[[i]], F)
  fit_filt$best_fit <-  best_fit
  p = plot_fit_density_mod(fit_filt, highlights = F, position = 'stack')
  p = cowplot::plot_grid(plotlist=p)
  
  title <- cowplot::ggdraw() + cowplot::draw_label(paste0("Likelihood: ", fit_filt$runs[[i]]$ICs$NLL), fontface='bold')
  p = cowplot::plot_grid(title, p, ncol=1, rel_heights=c(0.05, 1))
  
  ggsave(paste0(out_folder, 'k_', k, '_lambda_', lambda, "_summary_densities.png"),p,
         width = 25, height=25, dpi=150)
}


cell_names_atac <- Rcongas:::get_data(fit_filt) %>% filter(modality == "ATAC") %>%  pull(cell)  %>%  unique() %>%  sort()
cell_names_rna <- Rcongas:::get_data(fit_filt) %>% filter(modality == "RNA") %>%  pull(cell)  %>%  unique() %>%  sort()
gt_clusters_atac = celltypes %>% filter(modality == 'ATAC') %>% select(cell, type)
rownames(gt_clusters_atac) = gt_clusters_atac$cell
gt_clusters_atac$cell = NULL

gt_clusters_rna = celltypes %>% filter(modality == 'RNA') %>% select(cell, type)
rownames(gt_clusters_rna) = gt_clusters_rna$cell
gt_clusters_rna$cell = NULL

tmp = sapply(fit_filt$runs, function(w) {
  
  inferred_clusters_atac =  w$inferred_params$assignment_atac
  names(inferred_clusters_atac) = cell_names_atac
  inferred_clusters_atac = data.frame(cluster = inferred_clusters_atac[rownames(gt_clusters_atac)])
  
  inferred_clusters_rna =  w$inferred_params$assignment_rna
  names(inferred_clusters_rna) = cell_names_rna
  inferred_clusters_rna = data.frame(cluster = inferred_clusters_rna[rownames(gt_clusters_rna)])
  
  nmi_atac = aricode::NMI(gt_clusters_atac[,1], inferred_clusters_atac[,1])
  nmi_rna = aricode::NMI(gt_clusters_rna[,1], inferred_clusters_rna[,1])
  ari_atac = aricode::ARI(gt_clusters_atac[,1], inferred_clusters_atac[,1])
  ari_rna = aricode::ARI(gt_clusters_rna[,1], inferred_clusters_rna[,1])
  
  all_clusters_inferred = rbind(inferred_clusters_atac, inferred_clusters_rna)
  all_clusters_gt = rbind(gt_clusters_atac, gt_clusters_rna)
  
  nmi_total = aricode::NMI(all_clusters_gt[,1], all_clusters_inferred[,1])
  ari_total = aricode::ARI(all_clusters_gt[,1], all_clusters_inferred[,1])
  
  all_clusters_inferred = rbind(inferred_clusters_atac, inferred_clusters_rna)
  all_clusters_gt = rbind(gt_clusters_atac, gt_clusters_rna)
  
  nmi_total = aricode::NMI(all_clusters_gt[,1], all_clusters_inferred[,1])
  ari_total = aricode::ARI(all_clusters_gt[,1], all_clusters_inferred[,1])
  
  binded_groups_atac = cbind(inferred_clusters_atac %>% rename(inferred = cluster), gt_clusters_atac %>% rename(ground_truth = type)) %>% rownames_to_column('cell')
  binded_groups_atac = binded_groups_atac[gtools::mixedorder(binded_groups_atac$ground_truth),]
  # Add freq column to order the sankey plot
  binded_groups_atac$freq = seq(1, nrow(binded_groups_atac))
  library(ggalluvial)
  
  p = ggplot(data = binded_groups_atac,
             aes(axis1 = ground_truth, axis2 = inferred, y = freq)) +
    geom_alluvium(aes(fill = ground_truth)) +
    geom_stratum() +
    geom_text(stat = "stratum",
              aes(label = after_stat(stratum))) +
    # scale_x_discrete(limits = c("Survey", "Response"),
    # 				expand = c(0.15, 0.05)) +
    theme(axis.title=element_blank(), axis.text=element_blank(), axis.ticks=element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = "white"))
  
  k = w$hyperparameters$K
  lambda = w$hyperparameters$lambda
  #ggsave(p, filename = paste0(out_folder, 'ATAC_sankey_k_', k,  '_lambda_', lambda, '.png'), dpi = 150)
  
  return(list(nmi_rna = nmi_rna, nmi_atac = nmi_atac, nmi_total = nmi_total, ari_total = ari_total, ari_rna = ari_rna, ari_atac = ari_atac))
  
})

fit_filt$model_selection$nmi_rna = unlist(tmp[1,])
fit_filt$model_selection$nmi_atac = unlist(tmp[2,])
fit_filt$model_selection$nmi = unlist(tmp[3,])
fit_filt$model_selection$ari = unlist(tmp[4,])
fit_filt$model_selection$ari_rna = unlist(tmp[5,])
fit_filt$model_selection$ari_atac = unlist(tmp[6,])


scores = reshape2::melt(fit_filt$model_selection %>% 
                          select(NLL_rna, NLL_atac, NLL, BIC, K, lambda, nmi_atac, nmi_rna, ari_rna, ari_atac), id = c('K', 'lambda'))

scores$K = as.factor(scores$K)
IC_best = scores %>%
  group_by(variable) %>%
  filter(value == min(value)) 

p = scores %>% ggplot(aes(x = lambda, y = value, color = K)) +
  geom_point(
    data = IC_best,
    shape = 4,
    size = 4) + geom_line() +
  geom_point() +
  facet_wrap(~variable, scales = 'free') +
  #theme_linedraw() +
  scale_color_brewer(palette="Dark2") +
  theme_light() +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE, digits = 3))

ggsave(paste0(out_folder, 'NLL_BIC_ami_ari.pdf'), p)

# For every lambda I nees to select the best BIC and then consider the corresponding AMI value 
prova = fit_filt$model_selection %>% group_by(lambda, K) %>% arrange(BIC) %>% slice(1) %>% ungroup()
scores = reshape2::melt(prova %>% 
                          select(K, lambda, 
                                 nmi_atac, nmi_rna, 
                                 NLL_atac, NLL_rna
                                 #nmi
                          ), id = c('K', 'lambda'))

scores$K = as.factor(scores$K)
IC_best = bind_rows(scores %>% filter(variable %in% c('nmi_atac', 'nmi_rna', 'ari_rna', 'ari_atac')) %>%
                      group_by(variable) %>%
                      filter(value == max(value)), 
                    scores %>% filter(!(variable %in% c('nmi_atac', 'nmi_rna', 'ari_rna', 'ari_atac'))) %>%
                      group_by(variable) %>% filter(value == min(value)))

p = scores %>% ggplot(aes(x = lambda, y = value, color = K)) +
  geom_point(
    data = IC_best,
    shape = 4,
    size = 4) + 
  geom_line() +
  geom_point() +
  facet_wrap(~variable, scales = 'free') +
  #theme_linedraw() +
  scale_color_brewer(palette="Dark2") +
  theme_light() +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE, digits = 3))

ggsave(paste0(out_folder, 'ami_nll.pdf'), p)





