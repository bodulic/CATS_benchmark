#!/usr/bin/env Rscript
#Script for generating Figure 6 and Figure 7
#Loading the required packages
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggsci))
suppressPackageStartupMessages(library(cowplot))

#Importing merged controlled simulated assembly CATS-rb scores from file
controlled_simulated_assembly_scores <- fread("merged_cats_rb_controlled_simulated_assembly_scores_for_figures_6_7.tsv")

#Extracting species names
controlled_simulated_assembly_scores[, "species" := sub("_CATS.*", "", cats_rb_run)]
controlled_simulated_assembly_scores[, "species" := sub("_", ". ", species, fixed = T)]
controlled_simulated_assembly_scores[, "species" := paste0(toupper(substr(species, 1, 1)), substr(species, 2, nchar(species)))]

#Reshaping data to long format
controlled_simulated_assembly_scores_long <- melt(controlled_simulated_assembly_scores, measure.vars = colnames(controlled_simulated_assembly_scores)[c(-86, -87, -88)], variable.name = "mis_cov_assem", value.name = "score")

#Extracting coverage, mismatch rate, and assembler information
controlled_simulated_assembly_scores_long[, "coverage" := sub("^[^_]*_([^_]*_[^_]*)_.*", "\\1", mis_cov_assem)]
controlled_simulated_assembly_scores_long[, "coverage" := sub("_", "–", coverage, fixed = T)]
controlled_simulated_assembly_scores_long[coverage != "REF", "coverage" := paste0(coverage, "×")]

controlled_simulated_assembly_scores_long[, "mismatch" := sub("_.*", "", mis_cov_assem)]
controlled_simulated_assembly_scores_long[, "mismatch" := paste("ER", mismatch)]

controlled_simulated_assembly_scores_long[, "assembler" := sub(".*_", "", mis_cov_assem)]
controlled_simulated_assembly_scores_long[assembler == "RSP", "assembler" := "rnaSPAdes"]
controlled_simulated_assembly_scores_long[assembler == "TRI", "assembler" := "Trinity"]
controlled_simulated_assembly_scores_long[assembler == "IDB", "assembler" := "IDBA-tran"]
controlled_simulated_assembly_scores_long[assembler == "SOA", "assembler" := "SOAPdenovo-Trans"]

#Adjusting factor levels for species, coverage, mismatch, and assemblers
controlled_simulated_assembly_scores_long[, "species" := factor(species, levels = c("S. cerevisiae", "C. elegans", "D. melanogaster", "A. thaliana", "M. musculus", "H. sapiens"))]
controlled_simulated_assembly_scores_long[, "coverage" := factor(coverage, levels = c("1–4×", "5–10×", "11–20×", "21–30×", "31–40×", "41–50×", "51–100×"))]
controlled_simulated_assembly_scores_long[, "mismatch" := factor(mismatch, levels = c("ER 0.005", "ER 0.01", "ER 0.02"))]
controlled_simulated_assembly_scores_long[, "assembler" := factor(assembler, levels = c("rnaSPAdes", "Trinity", "IDBA-tran", "SOAPdenovo-Trans", "REF"))]

#Extracting CATS-rb runs on the full species datasets 
controlled_simulated_assembly_scores_full_dataset_long <- controlled_simulated_assembly_scores_long[grepl("sim$|sim1", cats_rb_run)]

#Plotting Figure 6a (distribution of CATS-rb relative transcript scores according to species and assembler) (controlled simulated assemblies)
fig_6a <- ggplot(controlled_simulated_assembly_scores_full_dataset_long[metric == "relative_transcript_score"], aes(x = assembler, y = score, color = assembler, fill = assembler)) +
 geom_boxplot(data = controlled_simulated_assembly_scores_full_dataset_long[metric == "relative_transcript_score" & assembler != "REF"], width = 0.8, size = 0.3, alpha = 0.75, color ="grey40", outlier.shape =  NA) +
 geom_point(position = position_jitter(width = 0.15, height = 0), size = 0.3, alpha = 0.8) +
 theme_minimal() +
 scale_fill_manual(values = c("#E64B35FF", "#4DBBD5FF", "#00A087FF","#3C5488FF", "#7F7F7F")) +
 scale_color_manual(values = c("#E64B35FF", "#4DBBD5FF", "#00A087FF","#3C5488FF", "#7F7F7F")) +
 theme(legend.position = "none")  +
 facet_wrap(. ~ species) +
 theme(strip.text = element_text(size = 7, face = "italic")) +
 theme(panel.grid = element_blank()) +
 theme(panel.background = element_rect(fill = NA, color = "grey70")) +
 xlab("Assembler") +
 ylab(expression(italic(R[t]))) +
 ylim(c(0, 1)) +
 theme(axis.title.x = element_text(size = 5.6))  +
 theme(axis.title.y = element_text(size = 7)) +
 theme(axis.text.x = element_text(size = 5.1, angle = 45, hjust = 0.9)) +
 theme(axis.text.y = element_text(size = 6.5)) 

ggsave(fig_6a, file = "Figure_6a.svg", height = 2.45, width = 3.3) 

#Plotting Supplementary Figure 4 (distribution of CATS-rb relative exon scores according to species and assembler) (controlled simulated assemblies)
supp_fig4 <- ggplot(controlled_simulated_assembly_scores_full_dataset_long[metric == "relative_exon_score"], aes(x = assembler, y = score, color = assembler, fill = assembler)) +
 geom_boxplot(data = controlled_simulated_assembly_scores_full_dataset_long[metric == "relative_exon_score" & assembler != "REF"], width = 0.8, size = 0.3, alpha = 0.75, color ="grey40", outlier.shape =  NA) +
 geom_point(position = position_jitter(width = 0.15, height = 0), size = 0.3, alpha = 0.8) +
 theme_minimal() +
 scale_fill_manual(values = c("#E64B35FF", "#4DBBD5FF", "#00A087FF","#3C5488FF", "#7F7F7F")) +
 scale_color_manual(values = c("#E64B35FF", "#4DBBD5FF", "#00A087FF","#3C5488FF", "#7F7F7F")) +
 theme(legend.position = "none")  +
 facet_wrap(. ~ species) +
 theme(strip.text = element_text(size = 5.5, face = "italic")) +
 theme(panel.grid = element_blank()) +
 theme(panel.background = element_rect(fill = NA, color = "grey70")) +
 xlab("Assembler") +
 ylab(expression(italic(R[e]))) +
 ylim(c(0, 1)) +
 theme(axis.title = element_text(size = 5.6)) +
 theme(axis.text.x = element_text(size = 4.5, angle = 40, hjust = 0.9)) +
 theme(axis.text.y = element_text(size = 4.5)) 

ggsave(supp_fig4, file = "Supplementary_Figure_4.tiff", height = 2.5, width = 4, dpi = 600, bg = "white") 

#Adding reference data for each error type (for Figure 6b amd Supplementary Figure 5)
controlled_simulated_assembly_scores_full_dataset_long_for_mis1 <- controlled_simulated_assembly_scores_full_dataset_long[assembler == "REF"]
controlled_simulated_assembly_scores_full_dataset_long_for_mis1[, "coverage" := "REF"]
controlled_simulated_assembly_scores_full_dataset_long_for_mis1[, "mismatch" := "ER 0.005"]

controlled_simulated_assembly_scores_full_dataset_long_for_mis2 <- controlled_simulated_assembly_scores_full_dataset_long[assembler == "REF"]
controlled_simulated_assembly_scores_full_dataset_long_for_mis2[, "coverage" := "REF"]
controlled_simulated_assembly_scores_full_dataset_long_for_mis2[, "mismatch" := "ER 0.01"]

controlled_simulated_assembly_scores_full_dataset_long_for_mis3 <- controlled_simulated_assembly_scores_full_dataset_long[assembler == "REF"]
controlled_simulated_assembly_scores_full_dataset_long_for_mis3[, "coverage" := "REF"]
controlled_simulated_assembly_scores_full_dataset_long_for_mis3[, "mismatch" := "ER 0.02"]

controlled_simulated_assembly_scores_full_dataset_long_for_fig6b_s5 <- rbindlist(list(controlled_simulated_assembly_scores_full_dataset_long[assembler != "REF"], controlled_simulated_assembly_scores_full_dataset_long_for_mis1, controlled_simulated_assembly_scores_full_dataset_long_for_mis2, controlled_simulated_assembly_scores_full_dataset_long_for_mis3))

#Plotting Figure 6b (distribution of CATS-rb relative transcript scores according to coverage and mismatch rate) (controlled simulated assemblies)
fig_6b <- ggplot(controlled_simulated_assembly_scores_full_dataset_long_for_fig6b_s5[metric == "relative_transcript_score"], aes(x = coverage, y = score, fill = coverage, color = coverage)) +
 geom_boxplot(width = 0.8, size = 0.3, alpha = 0.75, color ="grey40", outlier.shape =  NA) +
 geom_point(position = position_jitter(width = 0.15, height = 0), size = 0.15, alpha = 0.8) +
 theme_minimal() +
 scale_fill_manual(values = c("#36648B", "#4DBBD5FF", "#66C2A5", "#ABDDA4", "#F6C85F", "#F46D43", "#D53E4F", "#7F7F7F")) +
 scale_color_manual(values = c("#36648B", "#4DBBD5FF", "#66C2A5", "#ABDDA4", "#F6C85F", "#F46D43", "#D53E4F", "#7F7F7F")) +
 theme(legend.position = "none") +
 facet_grid(cols = vars(mismatch)) +
 theme(strip.text = element_text(size = 6)) +
 theme(panel.grid = element_blank()) +
 theme(panel.background = element_rect(fill = NA, color = "grey70")) +
 xlab("Coverage") +
 ylab(expression(italic(R[t]))) +
 ylim(c(0, 1)) +
 theme(axis.title.x = element_text(size = 6)) +
 theme(axis.title.y = element_text(size = 7)) +
 theme(axis.text.x = element_text(size = 5.6, angle = 45, hjust = 0.9)) +
 theme(axis.text.y = element_text(size = 6.5)) 

ggsave(fig_6b, file = "Figure_6b.svg", height = 1.4, width = 3.3) 

#Plotting Supplementary Figure 5 (distribution of CATS-rb relative exon scores according to mismatch rate and coverage) (controlled simulated assemblies)
supp_fig5 <- ggplot(controlled_simulated_assembly_scores_full_dataset_long_for_fig6b_s5[metric == "relative_exon_score"], aes(x = coverage, y = score, fill = coverage, color = coverage)) +
 geom_boxplot(width = 0.8, size = 0.3, alpha = 0.75, color ="grey40", outlier.shape =  NA) +
 geom_point(position = position_jitter(width = 0.15, height = 0), size = 0.15, alpha = 0.8) +
 theme_minimal() +
 scale_fill_manual(values = c("#36648B", "#4DBBD5FF", "#66C2A5", "#ABDDA4", "#F6C85F", "#F46D43", "#D53E4F", "#7F7F7F")) +
 scale_color_manual(values = c("#36648B", "#4DBBD5FF", "#66C2A5", "#ABDDA4", "#F6C85F", "#F46D43", "#D53E4F", "#7F7F7F")) +
 theme(legend.position = "none") +
 facet_grid(cols = vars(mismatch)) +
 theme(strip.text = element_text(size = 5.25)) +
 theme(panel.grid = element_blank()) +
 theme(panel.background = element_rect(fill = NA, color = "grey70")) +
 xlab("Coverage") +
 ylab(expression(italic(R[e]))) +
 ylim(c(0, 1)) +
 theme(axis.title = element_text(size = 5.5)) +
 theme(axis.text.x = element_text(size = 4.5, angle = 45, hjust = 0.9)) +
 theme(axis.text.y = element_text(size = 4.5)) 

ggsave(supp_fig5, file = "Supplementary_Figure_5.tiff", height = 1.5, width = 4, dpi = 600, bg = "white") 

#Plotting CATS-rb relative vs. annotation-based score scatterplots (function)
subset_sim_ids <- paste0("sim", 1 : 6)
plot_titles <- c("All assemblies", "51–100×", "31–40×", "11–20×", "5–10×", "1–4×")

plot_cats_rb_score_correlation <- function(rel_metric, annot_metric, output_file) {
 plot_list <- list()
 x_label <- if(grepl("transcript", annot_metric, fixed = T)) expression(A[t]) else expression(A[e])
 y_label <- if(grepl("transcript", rel_metric, fixed = T)) expression(R[t]) else expression(R[e])
 bg_color <- if (grepl("tiff", output_file)) "white" else "transparent"
 
 for (i in 1 : length(subset_sim_ids)) {
  sim_id <- subset_sim_ids[i]
  plot_title <- plot_titles[i]
  
  annot_score_table <- controlled_simulated_assembly_scores_long[metric == annot_metric, .(species, mis_cov_assem, score)]
  rel_score_table <- controlled_simulated_assembly_scores_long[grepl(sim_id, cats_rb_run, fixed = T) & metric == rel_metric, .(species, mis_cov_assem, score)]
  
  merged_score_table <- merge(rel_score_table, annot_score_table, by = c("species", "mis_cov_assem"))
  merged_score_table[, species := factor(species, levels = c("S. cerevisiae", "C. elegans", "D. melanogaster", "A. thaliana", "M. musculus", "H. sapiens"))]
    
  cor_coef <- sprintf("%.2f", cor(merged_score_table[, score.x], merged_score_table[, score.y], use = "pairwise.complete.obs", method = "spearman"))
  score_scatterplot <- ggplot(data = merged_score_table, aes(x = score.x, y = score.y, color = species)) +
   geom_point(size = 0.015, alpha = 0.8) +
   theme_minimal() +
   scale_color_manual(values = c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF", "#F6C85F")) +
   theme(legend.position = "none") +
   theme(panel.grid = element_blank()) +
   theme(panel.background = element_rect(fill = NA, color = "grey70")) +
   xlab(x_label) +
   ylab(y_label) +
   scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0, 1)) +
   scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0, 1)) +
   theme(axis.title = element_text(size = 7, face = "italic")) +
   theme(axis.text.x = element_text(size = 7, angle = 45, hjust = 0.9)) +
   theme(axis.text.y = element_text(size = 7)) +
   annotate("text", x = 0, y = 0.97, label = paste0("italic(\u03C1) == ", cor_coef), parse = T, size = 2.3, hjust = 0) +
   ggtitle(plot_title) +
   theme(plot.title = element_text(size = 7))
  
  plot_list[[i]] <- score_scatterplot
 }
 if(rel_metric == "relative_transcript_score") {
  final_plot <- plot_grid(plotlist = plot_list, nrow = 2)
 } else {
  legend_plot <- score_scatterplot + 
   theme(legend.position = "bottom") +
   theme(legend.title = element_blank()) +
   theme(legend.text = element_text(size = 4)) +
   theme(legend.key.size = unit(0.4, 'cm')) +
   theme(legend.key.width  = unit(1, "mm")) +
   theme(legend.key.height = unit(1, "mm"))
   
  legend <- get_legend(legend_plot)
  combined_plot <- plot_grid(plotlist = plot_list, nrow = 2)
  final_plot <- plot_grid(combined_plot, legend, ncol = 1, rel_heights = c(1, 0.12))
 }
 ggsave(final_plot, file = output_file, height = 3.2, width = 3.7, dpi = 600, bg = bg_color)
}

#Plotting Figure 6c (correlation of CATS-rb relative and annotation-based transcript scores in six analysed assembly subsets) (controlled simulated assemblies)
plot_cats_rb_score_correlation(rel_metric = "relative_transcript_score", annot_metric = "annotation_based_transcript_score", output_file = "Figure_6c.svg")

#Plotting Supplementary Figure 6 (correlation of CATS-rb relative and annotation-based exon scores in six analysed assembly subsets) (controlled simulated assemblies)
plot_cats_rb_score_correlation(rel_metric = "relative_exon_score", annot_metric = "annotation_based_exon_score", output_file = "Supplementary_Figure_6.tiff")

#Importing merged realistically simulated assembly CATS-rb scores from file
realistically_simulated_assembly_scores <- fread("merged_cats_rb_realistically_simulated_assembly_scores_for_figure6.tsv")

#Extracting species names
realistically_simulated_assembly_scores[, "species" := sub("_CATS.*", "", cats_rb_run)]
realistically_simulated_assembly_scores[, "species" := sub("_", ". ", species, fixed = T)]
realistically_simulated_assembly_scores[, "species" := paste0(toupper(substr(species, 1, 1)), substr(species, 2, nchar(species)))]

#Reshaping data to long format
realistically_simulated_assembly_scores_long <- melt(realistically_simulated_assembly_scores, measure.vars = colnames(realistically_simulated_assembly_scores)[c(-66, -67, -68)], variable.name = "cov_rep_assem", value.name = "score")

realistically_simulated_assembly_scores_long[, "coverage" := sub("_.*", "", cov_rep_assem)]
realistically_simulated_assembly_scores_long[coverage != "REF", "coverage" := paste0(coverage, "×")]

realistically_simulated_assembly_scores_long[, "assembler" := sub(".*_", "", cov_rep_assem)]
realistically_simulated_assembly_scores_long[assembler == "RSP", "assembler" := "rnaSPAdes"]
realistically_simulated_assembly_scores_long[assembler == "TRI", "assembler" := "Trinity"]
realistically_simulated_assembly_scores_long[assembler == "IDB", "assembler" := "IDBA-tran"]
realistically_simulated_assembly_scores_long[assembler == "SOA", "assembler" := "SOAPdenovo-Trans"]

#Adjusting factor levels for species, coverage, mismatch, and assemblers
realistically_simulated_assembly_scores_long[, "species" := factor(species, levels = c("S. cerevisiae", "C. elegans", "D. melanogaster", "A. thaliana", "M. musculus", "H. sapiens"))]
realistically_simulated_assembly_scores_long[, "coverage" := factor(coverage, levels = c("20×", "40×", "60×", "80×"))]
realistically_simulated_assembly_scores_long[, "assembler" := factor(assembler, levels = c("rnaSPAdes", "Trinity", "IDBA-tran", "SOAPdenovo-Trans", "REF"))]

realistically_simulated_assembly_scores_full_dataset_long_for_asse1 <- realistically_simulated_assembly_scores_long[assembler == "REF"]
realistically_simulated_assembly_scores_full_dataset_long_for_asse1[, "coverage" := "REF"]
realistically_simulated_assembly_scores_full_dataset_long_for_asse1[, "assembler" := "rnaSPAdes"]

realistically_simulated_assembly_scores_full_dataset_long_for_asse2 <- realistically_simulated_assembly_scores_long[assembler == "REF"]
realistically_simulated_assembly_scores_full_dataset_long_for_asse2[, "coverage" := "REF"]
realistically_simulated_assembly_scores_full_dataset_long_for_asse2[, "assembler" := "Trinity"]

realistically_simulated_assembly_scores_full_dataset_long_for_asse3 <- realistically_simulated_assembly_scores_long[assembler == "REF"]
realistically_simulated_assembly_scores_full_dataset_long_for_asse3[, "coverage" := "REF"]
realistically_simulated_assembly_scores_full_dataset_long_for_asse3[, "assembler" := "IDBA-tran"]

realistically_simulated_assembly_scores_full_dataset_long_for_asse4 <- realistically_simulated_assembly_scores_long[assembler == "REF"]
realistically_simulated_assembly_scores_full_dataset_long_for_asse4[, "coverage" := "REF"]
realistically_simulated_assembly_scores_full_dataset_long_for_asse4[, "assembler" := "SOAPdenovo-Trans"]

realistically_simulated_assembly_scores_full_dataset_long_for_fig6d_s5 <- rbindlist(list(realistically_simulated_assembly_scores_long[assembler != "REF"], realistically_simulated_assembly_scores_full_dataset_long_for_asse1, realistically_simulated_assembly_scores_full_dataset_long_for_asse2, realistically_simulated_assembly_scores_full_dataset_long_for_asse3, realistically_simulated_assembly_scores_full_dataset_long_for_asse4))

#Plotting Figure 6d (distribution of CATS-rb relative transcript scores according to assembler and coverage) (realistically simulated assemblies)
fig_6d <- ggplot(realistically_simulated_assembly_scores_full_dataset_long_for_fig6d_s5[metric == "relative_transcript_score"], aes(x = coverage, y = score, fill = coverage, color = coverage)) +
 geom_boxplot(width = 0.8, size = 0.3, alpha = 0.75, color ="grey40", outlier.shape =  NA) +
 geom_point(position = position_jitter(width = 0.15, height = 0), size = 0.15, alpha = 0.8) +
 theme_minimal() +
 scale_fill_manual(values = c( "#4DBBD5FF", "#ABDDA4", "#F6C85F", "#E64B35FF", "#7F7F7F")) +
 scale_color_manual(values = c( "#4DBBD5FF", "#ABDDA4", "#F6C85F", "#E64B35FF", "#7F7F7F")) +
 theme(legend.position = "none") +
 facet_grid(cols = vars(assembler)) +
 theme(strip.text = element_text(size = 5.18)) +
 theme(panel.grid = element_blank()) +
 theme(panel.background = element_rect(fill = NA, color = "grey70")) +
 xlab("Coverage") +
 ylab(expression(italic(R[t]))) +
 ylim(c(0, 1)) +
 theme(axis.text.x = element_text(size = 6.25, angle = 45, hjust = 0.9)) +
 theme(axis.text.y = element_text(size = 6.5)) +
 theme(axis.title.x = element_text(size = 6.18)) +
 theme(axis.title.y = element_text(size = 7))

ggsave(fig_6d, file = "Figure_6d.svg", height = 1.5, width = 4) 

#Plotting Supplementary Figure 7 (distribution of CATS-rb relative exon scores according to coverage and mismatch rate) (realistically simulated assemblies)
supp_fig7 <- ggplot(realistically_simulated_assembly_scores_full_dataset_long_for_fig6d_s5[metric == "relative_exon_score"], aes(x = coverage, y = score, fill = coverage, color = coverage)) +
 geom_boxplot(width = 0.8, size = 0.3, alpha = 0.75, color ="grey40", outlier.shape =  NA) +
 geom_point(position = position_jitter(width = 0.15, height = 0), size = 0.15, alpha = 0.8) +
 theme_minimal() +
 scale_fill_manual(values = c( "#4DBBD5FF", "#ABDDA4", "#F6C85F", "#E64B35FF", "#7F7F7F")) +
 scale_color_manual(values = c( "#4DBBD5FF", "#ABDDA4", "#F6C85F", "#E64B35FF", "#7F7F7F")) +
 theme(legend.position = "none") +
 facet_grid(cols = vars(assembler)) +
 theme(strip.text = element_text(size = 5.25)) +
 theme(panel.grid = element_blank()) +
 theme(panel.background = element_rect(fill = NA, color = "grey70")) +
 xlab("Coverage") +
 ylab(expression(italic(R[e]))) +
 ylim(c(0, 1)) +
 theme(axis.title = element_text(size = 5.5))
 theme(axis.text = element_text(size = 4.5)) +

ggsave(supp_fig7, file = "Supplementary_Figure_7.tiff", height = 1.5, width = 4.4, dpi = 600, bg = "white") 

#Importing merged realistically simulated assembly CATS-rb scores from file (per-library runs)
realistically_simulated_assembly_scores_per_library <- fread("merged_cats_rb_realistically_simulated_assembly_scores_per_library_for_figure6.tsv")

#Extracting species names
realistically_simulated_assembly_scores_per_library[, "species" := sub("^([^_]*_[^_]*).*", "\\1", cats_rb_run)]
realistically_simulated_assembly_scores_per_library[, "species" := sub("_", ". ", species, fixed = T)]
realistically_simulated_assembly_scores_per_library[, "species" := paste0(toupper(substr(species, 1, 1)), substr(species, 2, nchar(species)))]

#Adjusting factor levels for species
realistically_simulated_assembly_scores_per_library[, "species" := factor(species, levels = c("S. cerevisiae", "C. elegans", "D. melanogaster", "A. thaliana", "M. musculus", "H. sapiens"))]

#Reshaping data
realistically_simulated_assembly_scores_per_library_long <- melt(realistically_simulated_assembly_scores_per_library, measure.vars = colnames(realistically_simulated_assembly_scores_per_library)[1 : 4], variable.name = "assembler", value.name = "score")
realistically_simulated_assembly_scores_per_library_long_metric_wide <- dcast(realistically_simulated_assembly_scores_per_library_long, cats_rb_run + species + assembler ~ metric, value.var = "score")

#Plotting Figure 6e (Correlation of CATS-rb relative and annotation-based transcript/exon scores) (realistically simulated assemblies)
rel_annot_tr_score_cor_coeff <- sprintf("%.2f", cor(realistically_simulated_assembly_scores_per_library_long_metric_wide[, relative_transcript_score], realistically_simulated_assembly_scores_per_library_long_metric_wide[, annotation_based_transcript_score], method = "spearman"))
rel_annot_tr_score_scatterplot <- ggplot(data = realistically_simulated_assembly_scores_per_library_long_metric_wide, aes(x = annotation_based_transcript_score, y = relative_transcript_score, color = species)) +
 geom_point(size = 0.01, alpha = 0.8) +
 theme_minimal() +
 scale_color_manual(values = c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF", "#F6C85F")) +
 theme(legend.position = "none") +
 theme(panel.grid = element_blank()) +
 theme(panel.background = element_rect(fill = NA, color = "grey70")) +
 xlab(expression(italic(A[t]))) +
 ylab(expression(italic(R[t]))) +
 scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0, 1)) +
 scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0, 1)) +
 theme(axis.title = element_text(size = 7, face = "italic")) +
 theme(axis.text = element_text(size = 7)) +
 annotate("text", x = 0, y = 0.97, label = paste0("italic(\u03C1) == ", rel_annot_tr_score_cor_coeff), parse = T, size = 2.3, hjust = 0)

rel_annot_ex_cor_score_coeff <- sprintf("%.2f", cor(realistically_simulated_assembly_scores_per_library_long_metric_wide[, relative_exon_score], realistically_simulated_assembly_scores_per_library_long_metric_wide[, annotation_based_exon_score], method = "spearman"))
rel_annot_ex_score_scatterplot <- ggplot(data = realistically_simulated_assembly_scores_per_library_long_metric_wide, aes(x = annotation_based_exon_score, y = relative_exon_score, color = species)) +
 geom_point(size = 0.01, alpha = 0.8) +
 theme_minimal() +
 scale_color_manual(values = c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF", "#F6C85F")) +
 theme(legend.position = "none") +
 theme(panel.grid = element_blank()) +
 theme(panel.background = element_rect(fill = NA, color = "grey70")) +
 xlab(expression(italic(A[e]))) +
 ylab(expression(italic(R[e]))) +
 scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0, 1)) +
 scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0, 1)) +
 theme(axis.title = element_text(size = 7, face = "italic")) +
 theme(axis.text = element_text(size = 7)) +
 annotate("text", x = 0, y = 0.97, label = paste0("italic(\u03C1) == ", rel_annot_ex_cor_score_coeff), parse = T, size = 2.3, hjust = 0)

fig_6e <- plot_grid(rel_annot_tr_score_scatterplot, rel_annot_ex_score_scatterplot, nrow = 1)
ggsave(fig_6e, file = "Figure_6e.svg", height = 1.5, width = 3.8)

#Importing merged public assembly CATS-rb scores from file
simulated_public_scores <- fread("merged_cats_rb_public_assembly_scores_for_figure6.tsv")

#Extracting species names
simulated_public_scores[, "species" := sub("_SRR.*", "", cats_rb_run)]
simulated_public_scores[, "species" := sub("_", ". ", species, fixed = T)]
simulated_public_scores[, "species" := paste0(toupper(substr(species, 1, 1)), substr(species, 2, nchar(species)))]

#Adjusting factor levels for species
simulated_public_scores[, "species" := factor(species, levels = c("S. cerevisiae", "C. elegans", "D. melanogaster", "A. thaliana", "M. musculus", "H. sapiens"))]

#Reshaping data
simulated_public_scores_long <- melt(simulated_public_scores, measure.vars = colnames(simulated_public_scores)[1 : 4], variable.name = "assembler", value.name = "score")
simulated_public_scores_long_metric_wide <- dcast(simulated_public_scores_long, cats_rb_run + species + assembler ~ metric, value.var = "score")

#Plotting Figure 6f (Correlation of CATS-rb relative and annotation-based transcript/exon scores) (public assemblies)
rel_annot_tr_score_cor_coeff <- sprintf("%.2f", cor(simulated_public_scores_long_metric_wide[, relative_transcript_score], simulated_public_scores_long_metric_wide[, annotation_based_transcript_score], method = "spearman"))
rel_annot_tr_score_scatterplot <- ggplot(data = simulated_public_scores_long_metric_wide, aes(x = annotation_based_transcript_score, y = relative_transcript_score, color = species)) +
 geom_point(size = 0.02, alpha = 0.8) +
 theme_minimal() +
 scale_color_manual(values = c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF", "#F6C85F")) +
 theme(legend.position = "none") +
 theme(panel.grid = element_blank()) +
 theme(panel.background = element_rect(fill = NA, color = "grey70")) +
 xlab(expression(italic(A[t]))) +
 ylab(expression(italic(R[t]))) +
 scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0, 1)) +
 scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0, 1)) +
 theme(axis.title = element_text(size = 7, face = "italic")) +
 theme(axis.text = element_text(size = 7)) +
 annotate("text", x = 0, y = 0.97, label = paste0("italic(\u03C1) == ", rel_annot_tr_score_cor_coeff), parse = T, size = 2.3, hjust = 0)

rel_annot_ex_cor_score_coeff <- sprintf("%.2f", cor(simulated_public_scores_long_metric_wide[, relative_exon_score], simulated_public_scores_long_metric_wide[, annotation_based_exon_score], method = "spearman"))
rel_annot_ex_score_scatterplot <- ggplot(data = simulated_public_scores_long_metric_wide, aes(x = annotation_based_exon_score, y = relative_exon_score, color = species)) +
 geom_point(size = 0.02, alpha = 0.8) +
 theme_minimal() +
 scale_color_manual(values = c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF", "#F6C85F")) +
 theme(legend.position = "none") +
 theme(panel.grid = element_blank()) +
 theme(panel.background = element_rect(fill = NA, color = "grey70")) +
 xlab(expression(italic(A[e]))) +
 ylab(expression(italic(R[e]))) +
 scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0, 1)) +
 scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0, 1)) +
 theme(axis.title = element_text(size = 7, face = "italic")) +
 theme(axis.text = element_text(size = 7)) +
 annotate("text", x = 0, y = 0.97, label = paste0("italic(\u03C1) == ", rel_annot_ex_cor_score_coeff), parse = T, size = 2.3, hjust = 0)

fig_6f <- plot_grid(rel_annot_tr_score_scatterplot, rel_annot_ex_score_scatterplot, nrow = 1)
ggsave(fig_6f, file = "Figure_6f.svg", height = 1.5, width = 3.8)

#Importing merged CATS-rb simulated mutated assembly scores to file
simulated_mutated_assembly_scores <- fread("merged_cats_rb_simulated_mutated_assembly_scores_for_figure6.tsv")

#Reshaping data to long format
simulated_mutated_assembly_scores_long <- melt(simulated_mutated_assembly_scores, measure.vars = colnames(simulated_mutated_assembly_scores)[c(-39, -40)], variable.name = "assem_seed_level", value.name = "score")

#Denoting mutation level
simulated_mutated_assembly_scores_long[, "mut_level" := sub(".*_", "", assem_seed_level)]
simulated_mutated_assembly_scores_long[, "mut_level" := sub("0\\.", "", mut_level)]
simulated_mutated_assembly_scores_long[mut_level == "30", "mut_level" := "0"]

#Adjusting factor levels for mutation level
simulated_mutated_assembly_scores_long[, "mut_level" := factor(mut_level, levels = c("0", "1", "2", "3", "4", "5", "6"))]

#Plotting Figure 6g (distribution of CATS-rb relative transcript scores according to mutation level) (simulated multiplicative-mutation assemblies)
cats_rb_rel_tr_score_mut_boxplot <- ggplot(simulated_mutated_assembly_scores_long[metric == "relative_transcript_score"], aes(x = mut_level, y = score, color = mut_level, fill = mut_level)) +
 geom_boxplot(width = 0.8, size = 0.3, alpha = 0.75, color ="grey40", outlier.shape =  NA) +
 geom_point(position = position_jitter(width = 0.15, height = 0), size = 0.15, alpha = 0.8) +
 theme_minimal() +
 scale_fill_manual(values = c("#3C5488", "#4DBBD5", "#00A087", "#B09C85", "#F6C85F", "#F46D43", "#E64B35")) +
 scale_color_manual(values = c("#3C5488", "#4DBBD5", "#00A087", "#B09C85", "#F6C85F", "#F46D43", "#E64B35")) +
 theme(legend.position = "none")  +
 theme(panel.grid = element_blank()) +
 theme(panel.background = element_rect(fill = NA, color = "grey70")) +
 xlab("Mutation level") +
 ylab(expression(italic(R[t]))) +
 ylim(c(0, 1)) +
 theme(axis.title.x = element_text(size = 7)) +
 theme(axis.text = element_text(size = 7)) 

cats_rb_rel_ex_score_mut_boxplot <- ggplot(simulated_mutated_assembly_scores_long[metric == "relative_exon_score"], aes(x = mut_level, y = score, color = mut_level, fill = mut_level)) +
 geom_boxplot(width = 0.8, size = 0.3, alpha = 0.75, color ="grey40", outlier.shape =  NA) +
 geom_point(position = position_jitter(width = 0.15, height = 0), size = 0.15, alpha = 0.8) +
 theme_minimal() +
 scale_fill_manual(values = c("#3C5488", "#4DBBD5", "#00A087", "#B09C85", "#F6C85F", "#F46D43", "#E64B35")) +
 scale_color_manual(values = c("#3C5488", "#4DBBD5", "#00A087", "#B09C85", "#F6C85F", "#F46D43", "#E64B35")) +
 theme(legend.position = "none")  +
 theme(panel.grid = element_blank()) +
 theme(panel.background = element_rect(fill = NA, color = "grey70")) +
 xlab("Mutation level") +
 ylab(expression(italic(R[e]))) +
 ylim(c(0, 1)) +
 theme(axis.title = element_text(size = 7)) +
 theme(axis.text = element_text(size = 7)) 
 
fig_6g <- plot_grid(cats_rb_rel_tr_score_mut_boxplot, cats_rb_rel_ex_score_mut_boxplot, nrow = 1)
ggsave(fig_6g, file = "Figure_6g.svg", height = 1.5, width = 3.8)

#Importing CATS-rf assembly scores and mean transcript F-scores from file
simulated_assembly_cats_rf_f_scores <- fread("controlled_simulated_assembly_cats_rf_f_scores_for_figure7.tsv")

#Merging tables
controlled_sim_assembly_cats_rb_rf_f_scores <- merge(controlled_simulated_assembly_scores_long, simulated_assembly_cats_rf_f_scores, by = c("species", "coverage", "mismatch", "assembler"))

#Adjusting factor levels for species
controlled_sim_assembly_cats_rb_rf_f_scores[, "species" := factor(species, levels = c("S. cerevisiae", "C. elegans", "D. melanogaster", "A. thaliana", "M. musculus", "H. sapiens"))]

#Defining a function for plotting CATS-rb relative score vs. mean transcript F-score / CATS-rf assembly score
plot_cats_rb_f_cats_rf_score_correlation <- function(metric_type, x_var, x_label, y_label, output_file) {
 plot_list <- list()
 bg_color <- if (grepl("tiff", output_file)) "white" else "transparent"
 
 for (i in 1 : length(subset_sim_ids)) {
  sim_id <- subset_sim_ids[i]
  plot_title <- plot_titles[i]
    
  score_table_filt <- controlled_sim_assembly_cats_rb_rf_f_scores[grepl(sim_id, cats_rb_run, fixed = T) & metric == metric_type]
    
  cor_coef <- sprintf("%.2f", cor(score_table_filt[, score], score_table_filt[, ..x_var], use = "pairwise.complete.obs", method = "spearman"))
  score_scatterplot <- ggplot(data = score_table_filt, aes_string(x = x_var, y = "score", color = "species")) +
   geom_point(size = 0.015, alpha = 0.8) +
   theme_minimal() +
   scale_color_manual(values = c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF", "#F6C85F")) +
   theme(legend.position = "none") +
   theme(panel.grid = element_blank()) +
   theme(panel.background = element_rect(fill = NA, color = "grey70")) +
   xlab(x_label) +
   ylab(y_label) +
   scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0, 1)) +
   scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0, 1)) +
   theme(axis.title = element_text(size = 7)) +
   theme(axis.text.x = element_text(size = 7, angle = 45, hjust = 0.9)) +
   theme(axis.text.y = element_text(size = 7)) +
   annotate("text", x = 0, y = 0.97, label = paste0("italic(\u03C1) == ", cor_coef), parse = T, size = 2, hjust = 0) +
   ggtitle(plot_title) +
   theme(plot.title = element_text(size = 7))
    
  plot_list[[i]] <- score_scatterplot
 }
 
 if(grepl("transcript", metric_type)) {
  final_plot <- plot_grid(plotlist = plot_list, nrow = 2)
 } else {
  legend_plot <- score_scatterplot + 
   theme(legend.position = "bottom") +
   theme(legend.title = element_blank()) +
   theme(legend.text = element_text(size = 4)) +
   theme(legend.key.size = unit(0.32, 'cm')) +
   theme(legend.key.width  = unit(1, "mm")) +
   theme(legend.key.height = unit(1, "mm"))
  
  legend <- get_legend(legend_plot)
  combined_plot <- plot_grid(plotlist = plot_list, nrow = 2)
  final_plot <- plot_grid(combined_plot, legend, ncol = 1, rel_heights = c(1, 0.12))
 }
 ggsave(final_plot, file = output_file, height = 3, width = 4, dpi = 600, bg = bg_color)
}

#Plotting Figure 7a (Correlation of CATS-rb relative transcript scores with mean trasncripts F-scores in six analysed assembly subsets) (controlled simulated assemblies)
plot_cats_rb_f_cats_rf_score_correlation(metric_type = "relative_transcript_score", x_var = "transcript_f_score_mean", x_label = "Mean F-score", y_label = expression(R[t]), output_file = "Figure_7a.svg")

#Plotting Supplementary Figure 8 (Correlation of CATS-rb relative exon scores with mean trasncripts F-scores in six analysed assembly subsets) (controlled simulated assemblies)
plot_cats_rb_f_cats_rf_score_correlation(metric_type = "relative_exon_score", x_var = "transcript_f_score_mean", x_label = "Mean F-score", y_label = expression(R[e]), output_file = "Supplementary_figure_8.tiff")

#Plotting Figure 7b (Correlation of CATS-rb relative transcript scores with CATS-rf assembly scores in six analysed assembly subsets) (controlled simulated assemblies)
plot_cats_rb_f_cats_rf_score_correlation(metric_type = "relative_transcript_score", x_var = "cats_rf_assembly_score", x_label = expression(italic(S)), y_label = expression(R[t]), output_file = "Figure_7b.svg")

#Plotting Supplementary Figure 9 (Correlation of CATS-rb relative exon scores with CATS-rf assembly scores in six analysed assembly subsets) (controlled simulated assemblies)
plot_cats_rb_f_cats_rf_score_correlation(metric_type = "relative_exon_score", x_var = "cats_rf_assembly_score", x_label = expression(italic(S)), y_label = expression(R[e]), output_file = "Supplementary_figure_9.tiff")