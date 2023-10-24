# install.packages('~/OneDrive - University College London/dottorato/congas_git/rcongas/', repo = NULL, type = 'source')

machine = 'mac'
library(Rcongas)
library(stringi)
library(dplyr)
library(ggplot2)
library(stringr)
library(doParallel)
library(tidyverse)
if (machine == 'guascone') {
  dyn.load('/mnt/storage/dati_lucreziap/conda_envs/cn_atac/lib/libiconv.so.2')
  script_folder = '/mnt/storage/dati_lucreziap/congasp_analysis/'
  data_folder = '/mnt/storage/dati_lucreziap/Dropbox/CONGASp_data/'
} else if (machine == 'giacinto') {
  library(reticulate)
  reticulate::use_condaenv('/dati-raid/BimiB/share/conda_envs/congasp_reticulate/')
  script_folder = '/dati-raid/BimiB/lpatruno/ATAC/congasp_analysis'
  data_folder = '/dati-raid/BimiB/lpatruno/ATAC/CONGASp_data/'
} else if (machine == 'docker') {
  library(reticulate)
  reticulate::use_condaenv('/home/root/anaconda3')
  script_folder = '/app/congasp_analysis'
  data_folder = '/app/CONGASp_data'
} else if (machine == 'mac'){
  script_folder = '/Users/lucrezia/Library/CloudStorage/OneDrive-UniversityCollegeLondon/dottorato/congasp_analysis/'
  data_folder = '/Users/lucrezia/Dropbox/CONGASp_data/'
}

source(paste0(script_folder, "/congas_plot_mod.R"))

setwd(paste0(data_folder, "/10x_multiome"))
sample = 'lymphoma'
lik_rna = 'NB'
lik_atac = 'NB'
prefix = paste0(sample, "/bimodal_rna_", lik_rna, '_atac_', lik_atac)
library(tidyr)
cellularity = NULL#0.7

#x_filt = readRDS(paste0(prefix, 'congas_data/3.rcongas_top30.rds'))
lambda = c(0.2,0.5,0.8)
lambda_try = lambda
k = 1:3

binom_limits = c(15,1000)
model = "BIC"
steps_try = c(10000)

lr = 0.001
temperature = 20#10

steps = steps_try[1]
lambda = seq(0.05, 0.95, 0.1)

out_folder_prefix = '/lowLR_10ksteps_lowDisp_highPropb_1ksteps_lowlambda/'
out_folder_data = paste0(prefix, "/congas_data/", out_folder_prefix)
out_folder = paste0(prefix, "/congas_figures/", out_folder_prefix)

if (!dir.exists(out_folder)) {dir.create(out_folder, recursive = T)}
if (!dir.exists(out_folder_data)) {dir.create(out_folder_data, recursive = T)}

# if (!file.exists(paste0(prefix, 'congas_data/joint_assignments/single_segments.rds'))) {
if (!file.exists(paste0(prefix, "/congas_data/lowLR_10ksteps_lowDisp/single_segments.rds"))) {
    obj = readRDS(paste0(prefix, 
        '/congas_data/4.rcongas_noSameOutliers.rds'))
    x_filt = Rcongas::segments_selector_congas(obj, binom_limits = c(15,1000))
    saveRDS(x_filt, 
        paste0(prefix, "/congas_data//lowLR_10ksteps_lowDisp/single_segments.rds"))
} else {
    x_filt = readRDS(paste0(prefix, "/congas_data//lowLR_10ksteps_lowDisp/single_segments.rds"))
}

# n = readRDS('lymphoma/bimodal_rna_NB_atac_NB/congas_data/nic/fit.rds')
# x_filt = Rcongas:::select_segments(obj, segment_ids = n$input$segmentation$segment_id)


new_prior = F
cellularity = NULL

x_filt = Rcongas::filter_segments(x_filt, ATAC_peaks = 1000)

hyperparams_filt <- auto_config_run(x_filt, k, 
                                    prior_cn=c(0.2, 0.6, 0.1, 0.05, 0.05),
                                    multiome = TRUE,
                                    purity = cellularity, init_importance = 0.6,
                                    normal_cells = T)


hyperparams_filt$binom_prior_limits = binom_limits

fit_filt <- Rcongas:::fit_congas(x_filt,
                                 K = k, 
                                 lambdas = 0.2,#lambda, 
                                 learning_rate = lr, 
                                 steps = 5000,
                                 model_parameters = hyperparams_filt, 
                                 model_selection = model,
                                 latent_variables = "G",
                                 temperature = temperature, 
                                 same_mixing = TRUE, 
                                 threshold = 0.0001)


# prova = Rcongas:::fit_congas_single_run(x_filt, hyperparams_filt, l = 0.5, K = 3, lr,
#                                 latent_variables = 'G',  steps = 5, temperature = 20, patience = 5,
#                                 threshold = 0.1, same_mixing = FALSE, CUDA=F)

saveRDS(fit_filt, file=paste0(out_folder_data, "fit.rds"))
#fit_filt = readRDS(paste0(out_folder_data, 'fit.rds'))
#pdf(paste0(out_folder, "density.pdf"))
pdf(paste0(out_folder, "summary.pdf"))
p = plot_fit(fit_filt, highlights = T, what = 'density')
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
p = plot_fit_density_mod(fit_filt, highlights = F, position = 'stack')
print(p)
dev.off()
