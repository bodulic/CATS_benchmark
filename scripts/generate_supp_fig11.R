#!/usr/bin/env Rscript
#Script for generating Supplementary Figure 11
#Loading the required packages
library(data.table)
library(ggplot2)
library(cowplot)

#Importing CATS-rf transcript scores from files
cats_rf_tr_scores <- fread("cats_rf_transcript_scores_for_supp_fig11.tsv")
cats_rf_tr_scores[, "filename" := sub("_CATS_rf_param_res_transcript_scores.tsv", "", filename)]

#Extracting parameter values from file names 
all_parameter_values <- sub(".*_RSP_", "", cats_rf_tr_scores[, filename])
parameter_values <- gregexpr("([A-Za-z])([0-9]+(?:\\.[0-9]+)?)", all_parameter_values, perl = T)
parameter_values <- regmatches(all_parameter_values, parameter_values)
all_letters <- unique(unlist(lapply(parameter_values, function(param) {sub("([A-Za-z]).*", "\\1", param)})))
parameter_matrix <- matrix(NA_real_, nrow = cats_rf_tr_scores[, .N], ncol = length(all_letters), dimnames = list(NULL, all_letters))

for (i in 1 : length(parameter_values)) {
 parameter <- sub("([A-Za-z]).*", "\\1", parameter_values[[i]])
 parameter_value  <- as.numeric(sub("[A-Za-z]", "", parameter_values[[i]]))
 parameter_matrix[i, parameter] <- parameter_value
}
cats_rf_tr_scores  <- cbind(cats_rf_tr_scores, as.data.table(parameter_matrix))

#Importing transcript F-scores from files
tr_f_scores <- fread("transcript_f_scores_supp_fig11.tsv")
tr_f_scores[, "filename" := sub("_f_scores", "", filename)]
colnames(tr_f_scores)[1] <- "transcript"

#Merging CATS-rf transcript scores and F-scores
cats_rf_tr_scores[, "filename" := sub("(_RSP).*", "\\1", filename)]
cats_rf_tr_scores <- merge(cats_rf_tr_scores, tr_f_scores, by = c("transcript", "filename"))
cats_rf_tr_scores[, "cov_params" := paste(k, z, e, w, sep = "_")]
cats_rf_tr_scores[, "acc_params" := paste(K, Z, E, sep = "_")]
cats_rf_tr_scores[, "le_params" := paste(x, X, c, sep = "_")]
cats_rf_tr_scores[, "int_params" := paste(a, b, sep = "_")]

#Calculating and plotting the correlation between coverage score component and F-scores
cats_rf_tr_scores_cov <- cats_rf_tr_scores[, .("cor_coef" = cor(coverage_score_component, f_score, method = "spearman", use = "pairwise.complete.obs")), by = c("filename", "cov_params", "k")]
x_axis_face <-  unique(cats_rf_tr_scores_cov[order(k), cov_params])
x_axis_face <- fifelse(x_axis_face=="10_3_0.8_2", "bold", "plain")

cov_comp_plot <- ggplot(data = cats_rf_tr_scores_cov, aes(x = reorder(cov_params, k), y = cor_coef)) +
 geom_boxplot(width = 0.8, size = 0.2, alpha = 0.4, fill = "#E64B35", color ="grey40", outlier.shape = NA) +
 geom_jitter(size = 0.001, color = "#E64B35", alpha = 0.3) +
 theme_minimal() +
 theme(panel.grid = element_blank()) +
 theme(panel.background = element_rect(fill = NA, color = "grey70")) +
 xlab(expression("Parameters: "*italic(k)~" "*italic(z)*" "*italic(E)[c]*" "*italic(W))) +
 ylab("Correlation with F-score") +
 ylim(c(0, 1)) +
 theme(axis.title = element_text(size = 4.5)) +
 theme(axis.text.x = element_text(angle = 45, hjust = 0.9, size = 3, face = x_axis_face)) +
 theme(axis.text.y = element_text(size = 3.5)) +
 ggtitle(expression(italic(S[c]))) +
 theme(plot.title = element_text(size = 6)) 

#Calculating and plotting the correlation between accuracy score component and F-scores
cats_rf_tr_scores_acc <- cats_rf_tr_scores[, .("cor_coef" = cor(accuracy_score_component, f_score, method = "spearman", use = "pairwise.complete.obs")), by = c("filename", "acc_params", "K", "Z")]
x_axis_face <-  unique(cats_rf_tr_scores_acc[order(K, Z), acc_params])
x_axis_face <- fifelse(x_axis_face=="10_0.98_0.1", "bold", "plain")

acc_comp_plot <- ggplot(data = cats_rf_tr_scores_acc, aes(x = reorder(acc_params, K), y = cor_coef)) +
 geom_boxplot(width = 0.8, size = 0.2, alpha = 0.4, fill = "#F39B7F", color ="grey40", outlier.shape = NA) +
 geom_jitter(size = 0.001, color = "#F39B7F", alpha = 0.5) +
 theme_minimal() +
 theme(panel.grid = element_blank()) +
 theme(panel.background = element_rect(fill = NA, color = "grey70")) +
 xlab(expression("Parameters: "*italic(K)~" "*italic(Z)*" "*italic(E)[a])) +
 ylab("Correlation with F-score") +
 ylim(c(0, 1)) +
 theme(axis.title = element_text(size = 4.5)) +
 theme(axis.text.x = element_text(angle = 45, hjust = 0.9, size = 3, face = x_axis_face)) +
 theme(axis.text.y = element_text(size = 3.5)) +
 ggtitle(expression(italic(S[a]))) +
 theme(plot.title = element_text(size = 6)) 

#Calculating and plotting the correlation between local error score component and F-scores
cats_rf_tr_scores_loc <- cats_rf_tr_scores[, .("cor_coef" = cor(local_fidelity_score_component, f_score, method = "spearman", use = "pairwise.complete.obs")), by = c("filename", "le_params", "x")]
x_axis_face <-  unique(cats_rf_tr_scores_loc[order(x), le_params])
x_axis_face <- fifelse(x_axis_face=="8_10_5", "bold", "plain")

le_comp_plot <- ggplot(data = cats_rf_tr_scores_loc, aes(x = reorder(le_params, x), y = cor_coef)) +
 geom_boxplot(width = 0.8, size = 0.2, alpha = 0.4, fill = "#00A087FF", color ="grey40", outlier.shape = NA) +
 geom_jitter(size = 0.001, color = "#00A087FF", alpha = 0.5) +
 theme_minimal() +
 theme(panel.grid = element_blank()) +
 theme(panel.background = element_rect(fill = NA, color = "grey70")) +
 xlab(expression("Parameters: "*italic(M)[1]*" "*italic(M)[2]*" "*italic(C))) +
 ylab("Correlation with F-score") +
 ylim(c(0, 1)) +
 theme(axis.title = element_text(size = 4.5)) +
 theme(axis.text.x = element_text(angle = 45, hjust = 0.9, size = 3, face = x_axis_face)) +
 theme(axis.text.y = element_text(size = 3.5)) +
 ggtitle(expression(italic(S[l]))) +
 theme(plot.title = element_text(size = 6)) 

#Calculating and plotting the correlation between integrity score component and F-scores
cats_rf_tr_scores_int <- cats_rf_tr_scores[, .("cor_coef" = cor(integrity_score_component, f_score, method = "spearman", use = "pairwise.complete.obs")), by = c("filename", "int_params", "a")]
x_axis_face <-  unique(cats_rf_tr_scores_int[order(a), int_params])
x_axis_face <- fifelse(x_axis_face=="7_0.5", "bold", "plain")

int_comp_plot <- ggplot(data = cats_rf_tr_scores_int, aes(x = reorder(int_params, a), y = cor_coef)) +
 geom_boxplot(width = 0.8, size = 0.2, alpha = 0.4, fill = "#3C5488", color ="grey40", outlier.shape = NA) +
 geom_jitter(size = 0.005, color = "#3C5488", alpha = 0.5) +
 theme_minimal() +
 theme(panel.grid = element_blank()) +
 theme(panel.background = element_rect(fill = NA, color = "grey70")) +
 xlab(expression("Parameters: "*italic(alpha)*" "*italic(beta))) +
 ylab("Correlation with F-score") +
 ylim(c(0, 1)) +
 theme(axis.title = element_text(size = 4.5)) +
 theme(axis.text.x = element_text(angle = 45, hjust = 0.9, size = 3, face = x_axis_face)) +
 theme(axis.text.y = element_text(size = 3.5)) +
 ggtitle(expression(italic(S[i]))) +
 theme(plot.title = element_text(size = 6)) 

#Plotting Supplementary Figure 11 (Correlation between transcript score components and F-scores across different CATS-rf parameter combinations)
supp_fig11 <- plot_grid(cov_comp_plot, acc_comp_plot, le_comp_plot, int_comp_plot, nrow = 2, align = "hv")
ggsave(supp_fig11, file = "Supplementary_Figure_11.tiff", height = 3, width = 4, dpi = 600, bg = "white")