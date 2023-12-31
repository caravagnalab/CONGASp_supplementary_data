

library(Rcongas)
library(stringi)
library(dplyr)
library(ggplot2)
library(stringr)
library(doParallel)
library(tidyverse)


script_folder = './'
data_folder = 'data/'

source(paste0(script_folder, "/utils/congas_plot_mod.R"))
source(paste0(script_folder, "/utils/plot_cna_clusters.R"))


setwd(paste0(data_folder, "/SNU601/"))

lik_rna = 'NB'
lik_atac = 'NB'
prefix = paste0('bimodal_rna_', lik_rna, '_atac_', lik_atac, "/")
library(tidyr)
 

x_filt = readRDS(paste0(prefix, "congas_data/5.rcongas_noOutliers.rds"))

cellularity = NULL


k = c(1:8)

binom_limits = c(40,1000)
model = "BIC"
steps_try = c(15000)

lr = 0.01
temperature = 20

steps = steps_try[1]
lambda = 0.5


out_folder_prefix = paste0('congas_fit')
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


saveRDS(fit_filt, file=paste0(out_folder_data, ".rds"))


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
