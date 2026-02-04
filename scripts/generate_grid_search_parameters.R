#!/usr/bin/env Rscript
#Script for generating random grid search parameters for the analysis of CATS-rf parameter robustness
#Loading the required packages
suppressPackageStartupMessages(library(data.table))

#Coverage score component parameter value definition
k_values <- seq(from = 5, to = 30, by = 1)
z_values <- seq(from = 1.5, to = 5, by = 0.5)
e_values <- seq(from = 0.1, to = 1, by = 0.05)
w_values <- seq(from = 1, to = 5, by = 0.5)

#Accuracy score component parameter value definition
K_values <- seq(from = 5, to = 30, by = 1)
Z_values <- seq(from = 0.8, to = 1, by = 0.01)
E_values <- seq(from = 0.005, to = 0.5, by = 0.005)

#Local error score component parameter value definition
x_values <- seq(from = 2, to = 12, by = 0.5)
X_values <- seq(from = 4, to = 16, by = 0.5)
c_values <- seq(from = 1, to = 10, by = 1)

#Integrity score component parameter value definition
a_values <- seq(from = 2, to = 10, by = 0.5)
b_values <- seq(from = 0.5, to = 7, by = 0.5)

#Sampling coverage score component parameters
k_sample <- sample(k_values, size = 20, replace = T)
z_sample <- sample(z_values, size = 20, replace = T)
e_sample <- sample(e_values, size = 20, replace = T)
w_sample <- sample(w_values, size = 20, replace = T)

#Sampling accuracy score component parameters
K_sample <- sample(K_values, size = 20, replace = T)
Z_sample <- sample(Z_values, size = 20, replace = T)
E_sample <- sample(E_values, size = 20, replace = T)

#Sampling local error score component parameters
x_sample <- sample(x_values, size = 20, replace = T)
X_sample <- sample(X_values, size = 20, replace = T)
c_sample <- sample(c_values, size = 20, replace = T)

#Sampling integrity score component parameters
a_sample <- sample(a_values, size = 20, replace = T)
b_sample <- sample(b_values, size = 20, replace = T)

#Binding the parameters in a single data.table
parameter_dt <- as.data.table(cbind(k_sample, z_sample, e_sample, w_sample, K_sample, Z_sample, E_sample, x_sample, X_sample, c_sample, a_sample, b_sample))

#Adding default parameters
parameter_dt <- rbind(parameter_dt, list(10, 3, 0.5, 1.5, 10, 0.98, 0.1, 8, 10, 5, 7, 0.5))

#Constructing CATS-rf commands
fmt <- function(x) {
 if (is.integer(x)) return(as.character(x))
 format(x, scientific = F, trim = T)
}

parameter_dt[, "tag" := sprintf("k%s_z%s_e%s_w%s_K%s_Z%s_E%s_x%s_X%s_c%s_a%s_b%s", fmt(k_sample), fmt(z_sample), fmt(e_sample), fmt(w_sample), fmt(K_sample), fmt(Z_sample), fmt(E_sample), fmt(x_sample), fmt(X_sample), fmt(c_sample), fmt(a_sample), fmt(b_sample))]
parameter_dt[, "cmd" := sprintf(paste("CATS_rf -t 20", "-k %s -z %s -e %s -w %s", "-K %s -Z %s -E %s", "-x %s -X %s", "-c %s -a %s -b %s", "-o ${transcriptome}_%s_CATS_rf_param_res", "-D ${transcriptome}_%s_CATS_rf_param_dir", "../${transcriptome} ../sample_01_1.fasta ../sample_01_2.fasta"), fmt(k_sample), fmt(z_sample), fmt(e_sample), fmt(w_sample), fmt(K_sample), fmt(Z_sample), fmt(E_sample), fmt(x_sample), fmt(X_sample), fmt(c_sample), fmt(a_sample), fmt(b_sample), tag, tag)]

#Writing CATS-rf commands to file
write.table(parameter_dt[, .(cmd)], file = "CATS_rf_grid_search_commands", sep = "\n", row.names = F, col.names = F, quote = F)