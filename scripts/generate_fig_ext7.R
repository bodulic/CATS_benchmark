#!/usr/bin/env Rscript
#Script for generating Supplementary Figure 7
#Loading the required packages
library(data.table)
library(ggplot2)
library(ggsci)

#Importing chimeric reference transcript names from files
ref_tr_names <- fread("chimeric_ref_transcriptome_names_for_figure_ext7")
setnames(ref_tr_names, c("transcript", "transcriptome"))

#Importing CATS-rb structurally inconsistent transcripts
cats_rb_str_inc_tr <- fread("str_inconsistent_transcripts_for_figure_ext7")
setnames(cats_rb_str_inc_tr, c("transcript", "cats_rb_run"))

#Calculating the number of true-positive transcripts
tp_N <- cats_rb_str_inc_tr[grepl("_chimeric_mut$", transcript), .N, by = "cats_rb_run"]

#Calculating the number of false-negative transcripts
total_chim_N <- ref_tr_names[grepl("_chimeric_mut$", transcript), .N, by = "transcriptome"]
fn_N  <- total_chim_N[, N := N - tp_N[, N]]

#Calculating the number of false-positive transcripts
fp_N <- cats_rb_str_inc_tr[grepl("_chimeric_mut$", transcript) == F, .N, by = "cats_rb_run"]

#Calculating the number of true-negative transcripts
total_nonchim_N <- ref_tr_names[grepl("_chimeric_mut$", transcript) == F, .N, by = "transcriptome"]
tn_N <- total_nonchim_N[, N := N - fp_N[, N]]

#Merging the results
chim_detection_results <- cbind(tp_N, fn_N[, N], fp_N[, N], tn_N[, N])
setnames(chim_detection_results, c("transcriptome", "true_positives", "false_negatives", "false_positives", "true_negatives"))

#Calculating sensitivity and specificity
chim_detection_results[, "sensitivity" := true_positives / (true_positives + false_negatives)]
chim_detection_results[, "specificity" := true_negatives / (true_negatives + false_positives)]

#Reshaping data to long format
chim_detection_results_longer <- melt(chim_detection_results, measure.vars = c("sensitivity", "specificity"),  variable.name = "metric", value.name = "metric_value")

#Denoting species
chim_detection_results_longer[, "species" := sub("^(([^_]+)_([^_]+))_.*", "\\1", transcriptome)]
chim_detection_results_longer[, "species" := sub("_", ". ", species)]
chim_detection_results_longer[, "species" := paste0(toupper(substr(species, 1, 1)), substr(species, 2, nchar(species)))]

#Adjusting factor levels for species
chim_detection_results_longer[, "species" := factor(species, levels = c("S. cerevisiae", "C. elegans", "D. melanogaster", "A. thaliana", "M. musculus", "H. sapiens"))]

#Plotting Extended data Figure 7 (performance of CATS-rb in classifying chimeric sequences as structurally inconsistent)
fig_ext7 <- ggplot(data = chim_detection_results_longer, aes(x = species, y = metric_value, fill = metric)) +
 geom_bar(stat = "identity", position = "dodge", width = 0.7, size = 0.3, color = "grey35", alpha = 0.85) +
 theme_minimal() +
 scale_fill_npg(labels = c("Sensitivity", "Specificity")) +
 theme(legend.position = "right", legend.margin = margin(0, 0, 0, 0)) +
 theme(legend.title = element_blank()) +
 theme(legend.text = element_text(size = 4)) +
 theme(legend.key.size = unit(0.22, 'cm')) +
 theme(panel.grid = element_blank()) +
 theme(panel.background = element_rect(fill = NA, color = "grey70")) +
 scale_y_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1), limits = c(0, 1)) +
 theme(axis.title = element_blank()) +
 theme(axis.text = element_text(size = 4)) +
 theme(axis.text.x = element_text(face = "italic"))

ggsave(fig_ext7, file = "Extended_data_figure_7.tiff", height = 1.1, width = 3.2, dpi = 600, bg = "white")