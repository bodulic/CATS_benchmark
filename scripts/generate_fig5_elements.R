#!/usr/bin/env Rscript
#Script for generating Figure 5
#Loading the required packages
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggsci))
suppressPackageStartupMessages(library(cowplot))

#Importing merged simulated assembly CATS-rb scores from file
simulated_assembly_scores <- fread("merged_cats_rb_simulated_assembly_scores_for_figure5.tsv")

#Extracting species names
simulated_assembly_scores[, "species" := sub("_CATS.*", "", cats_rb_run)]
simulated_assembly_scores[, "species" := sub("_", ". ", species)]
simulated_assembly_scores[, "species" := paste0(toupper(substr(species, 1, 1)), substr(species, 2, nchar(species)))]

#Reshaping data to long format
simulated_assembly_scores_long <- melt(simulated_assembly_scores, measure.vars = colnames(simulated_assembly_scores)[c(-86, -87, -88)], variable.name = "mis_cov_assem", value.name = "score")

#Extracting coverage, mismatch rate, and assembler information
simulated_assembly_scores_long[, "coverage" := sub("^[^_]*_([^_]*_[^_]*)_.*", "\\1", mis_cov_assem)]
simulated_assembly_scores_long[, "coverage" := sub("_", "-", coverage)]
simulated_assembly_scores_long[, "coverage" := paste0(coverage, "x")]

simulated_assembly_scores_long[, "mismatch" := sub("_.*", "", mis_cov_assem)]
simulated_assembly_scores_long[, "mismatch" := paste("ER", mismatch)]

simulated_assembly_scores_long[, "assembler" := sub(".*_", "", mis_cov_assem)]
simulated_assembly_scores_long[assembler == "RSP", "assembler" := "rnaSPAdes"]
simulated_assembly_scores_long[assembler == "TRI", "assembler" := "Trinity"]
simulated_assembly_scores_long[assembler == "IDB", "assembler" := "IDBA-tran"]
simulated_assembly_scores_long[assembler == "SOA", "assembler" := "SOAPdenovo-Trans"]
simulated_assembly_scores_long[assembler == "ref", "assembler" := "REF"]

#Adjusting factor levels for species, coverage, mismatch, and assemblers
simulated_assembly_scores_long[, "species" := factor(species, levels = c("S. cerevisiae", "C. elegans", "D. melanogaster", "A. thaliana", "M. musculus", "H. sapiens"))]
simulated_assembly_scores_long[, "coverage" := factor(coverage, levels = c("1-4x", "5-10x", "11-20x", "21-30x", "31-40x", "41-50x", "51-100x"))]
simulated_assembly_scores_long[, "mismatch" := factor(mismatch, levels = c("ER 0.005", "ER 0.01", "ER 0.02"))]
simulated_assembly_scores_long[, "assembler" := factor(assembler, levels = c("rnaSPAdes", "Trinity", "IDBA-tran", "SOAPdenovo-Trans", "REF"))]

#Extracting CATS-rb runs on the full species datasets 
simulated_assembly_scores_full_dataset_long <- simulated_assembly_scores_long[grepl("sim$|sim1", cats_rb_run)]

#Plotting Figure 5a (distribution of CATS-rb relative transcript scores according to species and assembler) (simulated assemblies)
fig_5a <- ggplot(simulated_assembly_scores_full_dataset_long[metric == "relative_transcript_score"], aes(x = assembler, y = score, color = assembler, fill = assembler)) +
 geom_boxplot(width = 0.8, size = 0.3, alpha = 0.75, color ="grey40", outlier.shape =  NA) +
 geom_point(position = position_jitter(width = 0.15, height = 0), size = 0.15, alpha = 0.8) +
 theme_minimal() +
 scale_color_npg() +
 scale_fill_npg() +
 theme(legend.position = "none")  +
 facet_wrap(. ~ species) +
 theme(strip.text = element_text(size = 5.5, face = "italic")) +
 theme(panel.grid = element_blank()) +
 theme(panel.background = element_rect(fill = NA, color = "grey70")) +
 xlab("Assembler") +
 ylab(expression(R[t])) +
 ylim(c(0, 1)) +
 theme(axis.title = element_text(size = 5.6))  +
 theme(axis.title.y = element_text(face = "italic")) +
 theme(axis.text.x = element_text(size = 3.5, angle = 40, hjust = 0.9)) +
 theme(axis.text.y = element_text(size = 4.5)) 

ggsave(fig_5a, file = "Figure_5a.svg", height = 2.5, width = 4) 

#Plotting Extended data figure 2 (distribution of CATS-rb relative exon scores according to species and assembler) (simulated assemblies)
fig_ext2 <- ggplot(simulated_assembly_scores_full_dataset_long[metric == "relative_exon_score"], aes(x = assembler, y = score, color = assembler, fill = assembler)) +
 geom_boxplot(width = 0.8, size = 0.3, alpha = 0.75, color ="grey40", outlier.shape =  NA) +
 geom_point(position = position_jitter(width = 0.15, height = 0), size = 0.15, alpha = 0.8) +
 theme_minimal() +
 scale_color_npg() +
 scale_fill_npg() +
 theme(legend.position = "none")  +
 facet_wrap(. ~ species) +
 theme(strip.text = element_text(size = 5.5, face = "italic")) +
 theme(panel.grid = element_blank()) +
 theme(panel.background = element_rect(fill = NA, color = "grey70")) +
 xlab("Assembler") +
 ylab(expression(R[e])) +
 ylim(c(0, 1)) +
 theme(axis.title = element_text(size = 5.6))  +
 theme(axis.title.y = element_text(face = "italic")) +
 theme(axis.text.x = element_text(size = 3.5, angle = 40, hjust = 0.9)) +
 theme(axis.text.y = element_text(size = 4.5)) 

ggsave(fig_ext2, file = "Extended_data_figure_2.tiff", height = 2.5, width = 4, dpi = 600, bg = "white") 

#Adding reference data for each error type (for Figure 5b amd Supplementary Figure 3)
simulated_reference_assembly_scores_full_dataset_for_mis1 <- simulated_assembly_scores_full_dataset_long[assembler == "REF"]
simulated_reference_assembly_scores_full_dataset_for_mis1[, "coverage" := "REF"]
simulated_reference_assembly_scores_full_dataset_for_mis1[, "mismatch" := "ER 0.005"]

simulated_reference_assembly_scores_full_dataset_for_mis2 <- simulated_assembly_scores_full_dataset_long[assembler == "REF"]
simulated_reference_assembly_scores_full_dataset_for_mis2[, "coverage" := "REF"]
simulated_reference_assembly_scores_full_dataset_for_mis2[, "mismatch" := "ER 0.01"]

simulated_reference_assembly_scores_full_dataset_for_mis3 <- simulated_assembly_scores_full_dataset_long[assembler == "REF"]
simulated_reference_assembly_scores_full_dataset_for_mis3[, "coverage" := "REF"]
simulated_reference_assembly_scores_full_dataset_for_mis3[, "mismatch" := "ER 0.02"]

simulated_assembly_scores_full_dataset_long_for_fig5b_s3 <- rbindlist(list(simulated_assembly_scores_full_dataset_long[assembler != "REF"], simulated_reference_assembly_scores_full_dataset_for_mis1, simulated_reference_assembly_scores_full_dataset_for_mis2, simulated_reference_assembly_scores_full_dataset_for_mis3))

#Plotting Figure 5b (distribution of CATS-rb relative transcript scores according to coverage and mismatch rate) (simulated assemblies)
fig_5b <- ggplot(simulated_assembly_scores_full_dataset_long_for_fig5b_s3[metric == "relative_transcript_score"], aes(x = coverage, y = score, fill = coverage, color = coverage)) +
 geom_boxplot(width = 0.8, size = 0.3, alpha = 0.75, color ="grey40", outlier.shape =  NA) +
 geom_point(position = position_jitter(width = 0.15, height = 0), size = 0.15, alpha = 0.8) +
 theme_minimal() +
 scale_fill_manual(values = c("#36648B", "#3288BD", "#66C2A5", "#ABDDA4", "#FDAE61", "#F46D43", "#D53E4F", "#DC0000")) +
 scale_color_manual(values = c("#36648B", "#3288BD", "#66C2A5", "#ABDDA4", "#FDAE61", "#F46D43", "#D53E4F", "#DC0000")) +
 theme(legend.position = "none") +
 facet_grid(cols = vars(mismatch)) +
 theme(strip.text = element_text(size = 5.25)) +
 theme(panel.grid = element_blank()) +
 theme(panel.background = element_rect(fill = NA, color = "grey70")) +
 xlab("Coverage") +
 ylab(expression(R[t])) +
 ylim(c(0, 1)) +
 theme(axis.text.y = element_text(size = 4.5)) +
 theme(axis.text.x = element_text(size = 4.5, angle = 45, hjust = 0.9)) +
 theme(axis.title = element_text(size = 5.5)) +
 theme(axis.title.y = element_text(face = "italic"))

ggsave(fig_5b, file = "Figure_5b.svg", height = 1.5, width = 4) 

#Plotting Extended data Figure 3 (distribution of CATS-rb relative exon scores according to mismatch rate and coverage) (simulated assemblies)
fig_ext3 <- ggplot(simulated_assembly_scores_full_dataset_long_for_fig5b_s3[metric == "relative_exon_score"], aes(x = coverage, y = score, fill = coverage, color = coverage)) +
 geom_boxplot(width = 0.8, size = 0.3, alpha = 0.75, color ="grey40", outlier.shape =  NA) +
 geom_point(position = position_jitter(width = 0.15, height = 0), size = 0.15, alpha = 0.8) +
 theme_minimal() +
 scale_fill_manual(values = c("#36648B", "#3288BD", "#66C2A5", "#ABDDA4", "#FDAE61", "#F46D43", "#D53E4F", "#DC0000")) +
 scale_color_manual(values = c("#36648B", "#3288BD", "#66C2A5", "#ABDDA4", "#FDAE61", "#F46D43", "#D53E4F", "#DC0000")) +
 theme(legend.position = "none") +
 facet_grid(cols = vars(mismatch)) +
 theme(strip.text = element_text(size = 5.25)) +
 theme(panel.grid = element_blank()) +
 theme(panel.background = element_rect(fill = NA, color = "grey70")) +
 xlab("Coverage") +
 ylab(expression(R[e])) +
 ylim(c(0, 1)) +
 theme(axis.text.y = element_text(size = 4.5)) +
 theme(axis.text.x = element_text(size = 4.5, angle = 45, hjust = 0.9)) +
 theme(axis.title = element_text(size = 5.5)) +
 theme(axis.title.y = element_text(face = "italic"))

ggsave(fig_ext3, file = "Extended_data_figure_3.tiff", height = 1.5, width = 4, dpi = 600, bg = "white") 

#Plotting CATS-rb relative vs. annotation-based score scatterplots (function)
subset_sim_ids <- paste0("sim", 1 : 6)
plot_titles <- c("All assemblies", "51-100x", "31-40x", "11-20x", "5-10x", "1-4x")

plot_cats_rb_score_correlation <- function(rel_metric, annot_metric, output_file) {
 plot_list <- list()
 x_label <- if(grepl("transcript", annot_metric)) expression(A[t]) else expression(A[e])
 y_label <- if(grepl("transcript", rel_metric)) expression(R[t]) else expression(R[e])
 bg_color <- if (grepl("Extended", output_file)) "white" else "transparent"
 
 for (i in 1 : length(subset_sim_ids)) {
  sim_id <- subset_sim_ids[i]
  plot_title <- plot_titles[i]
  
  annot_score_table <- simulated_assembly_scores_long[metric == annot_metric, .(species, mis_cov_assem, score)]
  rel_score_table <- simulated_assembly_scores_long[grepl(sim_id, cats_rb_run) & metric == rel_metric, .(species, mis_cov_assem, score)]
  
  merged_score_table <- merge(rel_score_table, annot_score_table, by = c("species", "mis_cov_assem"))
  merged_score_table[, species := factor(species, levels = c("S. cerevisiae", "C. elegans", "D. melanogaster", "A. thaliana", "M. musculus", "H. sapiens"))]
    
  cor_coef <- sprintf("%.2f", cor(merged_score_table[, score.x], merged_score_table[, score.y], use = "pairwise.complete.obs", method = "spearman"))
  score_scatterplot <- ggplot(data = merged_score_table, aes(x = score.x, y = score.y, color = species)) +
   geom_point(size = 0.015, alpha = 0.8) +
   theme_minimal() +
   scale_color_npg() +
   theme(legend.position = "none") +
   theme(panel.grid = element_blank()) +
   theme(panel.background = element_rect(fill = NA, color = "grey70")) +
   xlab(x_label) +
   ylab(y_label) +
   scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0, 1)) +
   scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0, 1)) +
   theme(axis.title = element_text(size = 5.5, face = "italic")) +
   theme(axis.text = element_text(size = 5)) +
   annotate("text", x = 0, y = 1, label = paste("r =", cor_coef), size = 1.5, hjust = 0) +
   ggtitle(plot_title) +
   theme(plot.title = element_text(size = 5))
  
  plot_list[[i]] <- score_scatterplot
 }
 final_plot <- plot_grid(plotlist = plot_list, nrow = 2)
 ggsave(output_file, plot = final_plot, height = 3, width = 4, dpi = 600, bg = bg_color)
}

#Plotting Figure 5c (correlation of CATS-rb relative and annotation-based transcript scores in six analysed assembly subsets) (simulated assemblies)
plot_cats_rb_score_correlation(rel_metric = "relative_transcript_score", annot_metric = "annotation_based_transcript_score", output_file = "Figure_5c.svg")

#Plotting Extended data figure 4 (correlation of CATS-rb relative and annotation-based exon scores in six analysed assembly subsets) (simulated assemblies)
plot_cats_rb_score_correlation(rel_metric = "relative_exon_score", annot_metric = "annotation_based_exon_score", output_file = "Extended_data_figure_4.tiff")

#Importing merged public assembly CATS-rb scores from file
simulated_public_scores <- fread("merged_cats_rb_public_assembly_scores_for_figure5.tsv", select = 2 : 7)

#Extracting species names
simulated_public_scores[, "species" := sub("_SRR.*", "", cats_rb_run)]
simulated_public_scores[, "species" := sub("_", ". ", species)]
simulated_public_scores[, "species" := paste0(toupper(substr(species, 1, 1)), substr(species, 2, nchar(species)))]

#Adjusting factor levels for species
simulated_public_scores[, "species" := factor(species, levels = c("S. cerevisiae", "C. elegans", "D. melanogaster", "A. thaliana", "M. musculus", "H. sapiens"))]

#Reshaping data
simulated_public_scores_long <- melt(simulated_public_scores, measure.vars = colnames(simulated_public_scores)[1 : 4], variable.name = "assembler", value.name = "score")
simulated_public_scores_long_metric_wide <- dcast(simulated_public_scores_long, cats_rb_run + species + assembler ~ metric, value.var = "score")

#Plotting Figure 5d (Correlation of CATS-rb relative and annotation-based transcript/exon scores) (public assemblies)
rel_annot_tr_score_cor_coeff <- sprintf("%.2f", cor(simulated_public_scores_long_metric_wide[, relative_transcript_score], simulated_public_scores_long_metric_wide[, annotation_based_transcript_score], method = "spearman"))
rel_annot_tr_score_scatterplot <- ggplot(data = simulated_public_scores_long_metric_wide, aes(x = annotation_based_transcript_score, y = relative_transcript_score, color = species)) +
 geom_point(size = 0.02, alpha = 0.8) +
 theme_minimal() +
 scale_color_npg() +
 theme(legend.position = "none") +
 theme(panel.grid = element_blank()) +
 theme(panel.background = element_rect(fill = NA, color = "grey70")) +
 xlab(expression(A[t])) +
 ylab(expression(R[t])) +
 scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0, 1)) +
 scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0, 1)) +
 theme(axis.title = element_text(size = 5.5, face = "italic")) +
 theme(axis.text = element_text(size = 5)) +
 annotate("text", x = 0, y = 1, label = paste("r =", rel_annot_tr_score_cor_coeff), hjust = 0, size = 1.5) 

rel_annot_ex_cor_score_coeff <- sprintf("%.2f", cor(simulated_public_scores_long_metric_wide[, relative_exon_score], simulated_public_scores_long_metric_wide[, annotation_based_exon_score], method = "spearman"))
rel_annot_ex_score_scatterplot <- ggplot(data = simulated_public_scores_long_metric_wide, aes(x = annotation_based_exon_score, y = relative_exon_score, color = species)) +
 geom_point(size = 0.02, alpha = 0.8) +
 theme_minimal() +
 scale_color_npg() +
 theme(legend.position = "none") +
 theme(panel.grid = element_blank()) +
 theme(panel.background = element_rect(fill = NA, color = "grey70")) +
 xlab(expression(A[e])) +
 ylab(expression(R[e])) +
 scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0, 1)) +
 scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0, 1)) +
 theme(axis.title = element_text(size = 5.5, face = "italic")) +
 theme(axis.text = element_text(size = 5)) +
 annotate("text", x = 0, y = 1, label = paste("r =", rel_annot_ex_cor_score_coeff), hjust = 0, size = 1.5) 

fig_5d <- plot_grid(rel_annot_tr_score_scatterplot, rel_annot_ex_score_scatterplot, nrow = 1)
ggsave(fig_5d, file = "Figure_5d.svg", height = 1.5, width = 3.8)

#Importing CATS-rf assembly scores and mean transcript F-scores from file
simulated_assembly_cats_rf_f_scores <- fread("simulated_assembly_cats_rf_f_scores_for_figure5.tsv")

#Merging tables
sim_assembly_cats_rb_rf_f_scores <- merge(simulated_assembly_scores_long, simulated_assembly_cats_rf_f_scores, by = c("species", "coverage", "mismatch", "assembler"))

#Adjusting factor levels for species
sim_assembly_cats_rb_rf_f_scores[, "species" := factor(species, levels = c("S. cerevisiae", "C. elegans", "D. melanogaster", "A. thaliana", "M. musculus", "H. sapiens"))]

#Defining a function for plotting CATS-rb relative score vs. mean transcript F-score / CATS-rf assembly score
plot_cats_rb_f_cats_rf_score_correlation <- function(metric_type, x_var, x_label, y_label, output_file) {
 plot_list <- list()
 bg_color <- if (grepl("Extended", output_file)) "white" else "transparent"
 
 for (i in 1 : length(subset_sim_ids)) {
  sim_id <- subset_sim_ids[i]
  plot_title <- plot_titles[i]
    
  score_table_filt <- sim_assembly_cats_rb_rf_f_scores[grepl(sim_id, cats_rb_run) & metric == metric_type]
  
  cor_coef <- sprintf("%.2f", cor(score_table_filt[, score], score_table_filt[, ..x_var], use = "pairwise.complete.obs", method = "spearman"))
  
  score_scatterplot <- ggplot(data = score_table_filt, aes_string(x = x_var, y = "score", color = "species")) +
   geom_point(size = 0.015, alpha = 0.8) +
   theme_minimal() +
   scale_color_npg() +
   theme(legend.position = "none") +
   theme(panel.grid = element_blank()) +
   theme(panel.background = element_rect(fill = NA, color = "grey70")) +
   xlab(x_label) +
   ylab(y_label) +
   scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0, 1)) +
   scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0, 1)) +
   theme(axis.title = element_text(size = 5.5)) +
   theme(axis.title.y = element_text(face = "italic")) +
   theme(axis.text = element_text(size = 5)) +
   annotate("text", x = 0, y = 1, label = paste("r =", cor_coef), size = 1.5, hjust = 0) +
   ggtitle(plot_title) +
   theme(plot.title = element_text(size = 5))
    
  plot_list[[i]] <- score_scatterplot
 }
 final_plot <- plot_grid(plotlist = plot_list, nrow = 2)
 ggsave(final_plot, file = output_file, height = 3, width = 4, dpi = 600, bg = bg_color)
}

#Plotting Figure 5e (Correlation of CATS-rb relative transcript scores with mean trasncripts F-scores in six analysed assembly subsets) (simulated assemblies)
plot_cats_rb_f_cats_rf_score_correlation(metric_type = "relative_transcript_score", x_var = "transcript_f_score_mean", x_label = "Mean F-score", y_label = expression(R[t]), output_file = "Figure_5e.svg")

#Plotting Extended data figure 5 (Correlation of CATS-rb relative exon scores with mean trasncripts F-scores in six analysed assembly subsets) (simulated assemblies)
plot_cats_rb_f_cats_rf_score_correlation(metric_type = "relative_exon_score", x_var = "transcript_f_score_mean", x_label = "Mean F-score", y_label = expression(R[e]), output_file = "Extended_data_figure5.tiff")

#Plotting Figure 5f (Correlation of CATS-rb relative transcript scores with CATS-rf assembly scores in six analysed assembly subsets) (simulated assemblies)
plot_cats_rb_f_cats_rf_score_correlation(metric_type = "relative_transcript_score", x_var = "cats_rf_assembly_score", x_label = expression(italic(S)), y_label = expression(R[t]), output_file = "Figure_5f.svg")

#Plotting Extended data figure 6 (Correlation of CATS-rb relative exon scores with CATS-rf assembly scores in six analysed assembly subsets) (simulated assemblies)
plot_cats_rb_f_cats_rf_score_correlation(metric_type = "relative_exon_score", x_var = "cats_rf_assembly_score", x_label = expression(italic(S)), y_label = expression(R[e]), output_file = "Extended_data_figure6.tiff")

#Importing merged CATS-rb simulated mutated assembly scores to file
simulated_mutated_assembly_scores <- fread("merged_cats_rb_simulated_mutated_assembly_scores_for_figure5.tsv")

#Reshaping data to long format
simulated_mutated_assembly_scores_long <- melt(simulated_mutated_assembly_scores, measure.vars = colnames(simulated_mutated_assembly_scores)[c(-39, -40)], variable.name = "assem_seed_level", value.name = "score")

#Denoting mutation level
simulated_mutated_assembly_scores_long[, "mut_level" := sub(".*_", "", assem_seed_level)]
simulated_mutated_assembly_scores_long[, "mut_level" := sub("0\\.", "", mut_level)]
simulated_mutated_assembly_scores_long[mut_level == "30", "mut_level" := "0"]

#Adjusting factor levels for mutation level
simulated_mutated_assembly_scores_long[, "mut_level" := factor(mut_level, levels = c("0", "1", "2", "3", "4", "5", "6"))]

#Plotting Figure 5g (distribution of CATS-rb relative transcript scores according to mutation level) (simulated multiplicative-mutation assemblies)
cats_rb_rel_tr_score_mut_boxplot <- ggplot(simulated_mutated_assembly_scores_long[metric == "relative_transcript_score"], aes(x = mut_level, y = score, color = mut_level, fill = mut_level)) +
 geom_boxplot(width = 0.8, size = 0.3, alpha = 0.75, color ="grey40", outlier.shape =  NA) +
 geom_point(position = position_jitter(width = 0.15, height = 0), size = 0.15, alpha = 0.8) +
 theme_minimal() +
 scale_fill_manual(values = c("#3C5488", "#4DBBD5", "#00A087", "#B09C85", "#FDAE61", "#F46D43", "#E64B35")) +
 scale_color_manual(values = c("#3C5488", "#4DBBD5", "#00A087", "#B09C85", "#FDAE61", "#F46D43", "#E64B35")) +
 theme(legend.position = "none")  +
 theme(panel.grid = element_blank()) +
 theme(panel.background = element_rect(fill = NA, color = "grey70")) +
 xlab("Mutation level") +
 ylab(expression(R[t])) +
 ylim(c(0, 1)) +
 theme(axis.title = element_text(size = 6.3)) +
 theme(axis.title.y = element_text(face = "italic")) +
 theme(axis.text = element_text(size = 5.8)) 

cats_rb_rel_ex_score_mut_boxplot <- ggplot(simulated_mutated_assembly_scores_long[metric == "relative_exon_score"], aes(x = mut_level, y = score, color = mut_level, fill = mut_level)) +
 geom_boxplot(width = 0.8, size = 0.3, alpha = 0.75, color ="grey40", outlier.shape =  NA) +
 geom_point(position = position_jitter(width = 0.15, height = 0), size = 0.15, alpha = 0.8) +
 theme_minimal() +
 scale_fill_manual(values = c("#3C5488", "#4DBBD5", "#00A087", "#B09C85", "#FDAE61", "#F46D43", "#E64B35")) +
 scale_color_manual(values = c("#3C5488", "#4DBBD5", "#00A087", "#B09C85", "#FDAE61", "#F46D43", "#E64B35")) +
 theme(legend.position = "none")  +
 theme(panel.grid = element_blank()) +
 theme(panel.background = element_rect(fill = NA, color = "grey70")) +
 xlab("Mutation level") +
 ylab(expression(R[e])) +
 ylim(c(0, 1)) +
 theme(axis.title = element_text(size = 6.3)) +
 theme(axis.title.y = element_text(face = "italic")) +
 theme(axis.text = element_text(size = 5.8)) 

fig_5g <- plot_grid(cats_rb_rel_tr_score_mut_boxplot, cats_rb_rel_ex_score_mut_boxplot, nrow = 1)
ggsave(fig_5g, file = "Figure_5g.svg", height = 1.5, width = 3.8)