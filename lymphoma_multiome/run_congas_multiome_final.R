# install.packages('~/OneDrive - University College London/dottorato/congas_git/rcongas/', repo = NULL, type = 'source')

machine = 'mac'
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

setwd(paste0(data_folder, "/10x_multiome"))
sample = 'lymphoma'
lik_rna = 'NB'
lik_atac = 'NB'
prefix = paste0(sample, "/bimodal_rna_", lik_rna, '_atac_', lik_atac)
library(tidyr)
cellularity = NULL
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

out_folder_prefix = '/congas_fit'
out_folder_data = paste0(prefix, "/congas_data/", out_folder_prefix)
out_folder = paste0(prefix, "/congas_figures/", out_folder_prefix)

if (!dir.exists(out_folder)) {dir.create(out_folder, recursive = T)}
if (!dir.exists(out_folder_data)) {dir.create(out_folder_data, recursive = T)}

x_filt = Rcongas::segments_selector_congas(x)

cellularity = NULL


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


saveRDS(fit_filt, file=paste0(out_folder_data, ".rds"))
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
