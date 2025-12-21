#!/usr/bin/env Rscript
#Script for generating Figure 2
#Loading the required packages
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggsci))
suppressPackageStartupMessages(library(ggcorrplot))
suppressPackageStartupMessages(library(cowplot))

#Importing merged controlled simulated transcript scores from file
controlled_simulated_transcript_scores <- fread("merged_controlled_simulated_transcript_scores_for_figure2.tsv")

#Extracting species names, coverage, mismatch rate, and assembler information
controlled_simulated_transcript_scores[, "assembly" := sub("controlled_", "", assembly, fixed = T)]
controlled_simulated_transcript_scores[, "assembly" := sub("_1_12345", "", assembly, fixed = T)]

controlled_simulated_transcript_scores[, "species" := sub("^sim_([^_]+_[^_]+)_.*", "\\1",  assembly)]
controlled_simulated_transcript_scores[, "species" := sub("_", ". ", species, fixed = T)]
controlled_simulated_transcript_scores[, "species" := paste0(toupper(substr(species, 1, 1)), substr(species, 2, nchar(species)))]

controlled_simulated_transcript_scores[, "coverage" := sub("^.*\\.fa_100_([^_]+_[^_]+)_[^_]+_[^_]+$", "\\1", assembly)]
controlled_simulated_transcript_scores[, "coverage" := sub("_", "-", coverage, fixed = T)]
controlled_simulated_transcript_scores[, "coverage" := paste0(coverage, "x")]

controlled_simulated_transcript_scores[, "mismatch" := sub(".*_([0-9.]+)_[^_]+$", "\\1", assembly)]
controlled_simulated_transcript_scores[, "mismatch" := paste("ER", mismatch)]

controlled_simulated_transcript_scores[, "assembler" := sub(".*_", "", assembly)]
controlled_simulated_transcript_scores[assembler == "RSP", "assembler" := "rnaSPAdes"]
controlled_simulated_transcript_scores[assembler == "TRI", "assembler" := "Trinity"]
controlled_simulated_transcript_scores[assembler == "IDB", "assembler" := "IDBA-tran"]
controlled_simulated_transcript_scores[assembler == "SOA", "assembler" := "SOAPdenovo-Trans"]

#Calculating and writing mean transcript F-scores and CATS-rf assembly scores to file (for Figures 5g and Figure 5h)
controlled_simulated_assem_cats_rf_f_scores <- controlled_simulated_transcript_scores[, .("cats_rf_assembly_score" = mean(cats_rf_transcript_score), "transcript_f_score_mean" = mean(transcript_f_score)), by = c("species", "coverage", "mismatch", "assembler")]
write.table(controlled_simulated_assem_cats_rf_f_scores, file = "controlled_simulated_assembly_cats_rf_f_scores_for_figure5.tsv", sep = "\t", row.names = F, col.names = T, quote = F)

#Adjusting factor levels for species, coverage, mismatch, and assemblers
controlled_simulated_transcript_scores[, "species" := factor(species, levels = c("S. cerevisiae", "C. elegans", "D. melanogaster", "A. thaliana", "M. musculus", "H. sapiens"))]
controlled_simulated_transcript_scores[, "coverage" := factor(coverage, levels = c("1-4x", "5-10x", "11-20x", "21-30x", "31-40x", "41-50x", "51-100x"))]
controlled_simulated_transcript_scores[, "mismatch" := factor(mismatch, levels = c("ER 0.005", "ER 0.01", "ER 0.02"))]
controlled_simulated_transcript_scores[, "assembler" := factor(assembler, levels = c("rnaSPAdes", "Trinity", "IDBA-tran", "SOAPdenovo-Trans"))]

#Plotting Figure 2a (distribution of CATS-rf scores according to species and assembler) (controlled simulated assemblies)
fig_2a <- ggplot(data = controlled_simulated_transcript_scores, aes(x = cats_rf_transcript_score, group = assembler, color = assembler)) +
 geom_density(aes(y = ..density.. / max(..density..)), adjust = 3, linewidth = 0.4) +
 theme_minimal() +
 scale_color_npg(name = "Assembler") +
 theme(legend.position = "bottom", legend.justification = "center", legend.margin = margin(0, 0, 0, 0), legend.box.margin = margin(-6, -6, -3, -6)) +
 theme(legend.title = element_text(size = 5.2)) +
 theme(legend.text = element_text(size = 4.8)) +
 theme(legend.key.size = unit(0.3, 'cm')) + 
 facet_wrap(. ~ species) +
 theme(strip.text = element_text(size = 6.7, face = "italic"))  +
 theme(panel.grid = element_blank())  +
 theme(panel.background = element_rect(fill = NA, color = "grey70")) +
 xlab(expression(S[t])) +
 ylab("Relative density") +
 theme(axis.title = element_text(size = 5.5)) +
 theme(axis.title.x = element_text(face = "italic")) +
 theme(axis.text = element_text(size = 5.3)) 
  
ggsave(fig_2a, file = "Figure_2a.svg", height = 3, width = 4.3)

#Plotting Figure 2b (distribution of CATS-rf scores according to coverage, mismatch rate, and assembler) (controlled simulated assemblies)
fig_2b <- ggplot(data = controlled_simulated_transcript_scores, aes(x = cats_rf_transcript_score, group = coverage, color = coverage)) +
 geom_density(aes(y = ..density.. / max(..density..)), adjust = 3, linewidth = 0.2) +
 theme_minimal() +
 scale_color_manual(values = c("#3C5488", "#4DBBD5", "#00A087", "#B09C85", "#FDAE61", "#F46D43", "#E64B35"), name = "Coverage") +
 theme(legend.position = "bottom", legend.justification = "center", legend.margin = margin(0, 0, 0, 0), legend.box.margin = margin(-6, -6, -3, -6)) +
 guides(colour = guide_legend(nrow = 1)) +
 theme(legend.title = element_text(size = 5.4)) +
 theme(legend.text = element_text(size = 5)) +
 theme(legend.key.size = unit(0.28, 'cm')) + 
 facet_grid(assembler ~ mismatch) +
 theme(strip.text.x = element_text(size = 5)) + 
 theme(strip.text.y = element_text(size = 4.5)) + 
 theme(panel.grid = element_blank()) +
 theme(panel.background = element_rect(fill = NA, color = "grey70")) +
 xlab(expression(S[t])) +
 ylab("Relative density") +
 theme(axis.title = element_text(size = 5.6)) +
 theme(axis.title.x = element_text(face = "italic")) +
 theme(axis.text = element_text(size = 5)) 

ggsave(fig_2b, file = "Figure_2b.svg", height = 4, width = 4)
  
#Calculating pairwise correlation between CATS-rf score components
assemblies <- unique(controlled_simulated_transcript_scores[, assembly])
pairwise_cats_rf_comp_corr_matrix_list <- list()

for (i in 1 : length(assemblies)) {
 controlled_simulated_transcript_scores_assembly <- controlled_simulated_transcript_scores[assembly == assemblies[i]]
 pairwise_cats_rf_comp_corr_matrix_list[[i]] <- cor(controlled_simulated_transcript_scores_assembly[, .(cats_rf_coverage_component, cats_rf_accuracy_component, cats_rf_local_fidelity_component, cats_rf_integrity_component)], use = "pairwise.complete.obs", method = "spearman")
}
pairwise_cats_rf_comp_corr_matrix_array <- simplify2array(pairwise_cats_rf_comp_corr_matrix_list)
pairwise_cats_rf_comp_corr_matrix_median <- apply(pairwise_cats_rf_comp_corr_matrix_array, c(1, 2), median)

#Plotting Figure 2c (correlation matrix of CATS-rf score components) (controlled simulated assemblies)
fig_2c <- ggcorrplot(pairwise_cats_rf_comp_corr_matrix_median, method = "square", type = "full", ggtheme = ggplot2::theme_minimal(), colors = c("#E64B35", "white", "#4DBBD5"), show.legend = F, lab = T, show.diag = T, lab_size = 2.1, digits = 2, tl.cex = 8, tl.srt = 0) +
 scale_x_discrete(labels = c(expression(S[c]), expression(S[a]), expression(S[l]), expression(S[i]))) + 
 scale_y_discrete(labels = c(expression(S[c]), expression(S[a]), expression(S[l]), expression(S[i]))) +
 theme(axis.text = element_text(face = "italic", size = 13))

ggsave(fig_2c, file = "figure_2c.svg", height = 1.6, width = 1.6)

#Calculating correlations of CATS-rf, RSEM-EVAL, and TransRate transcript scores and score components with transcript F-scores
metrics <- c("cats_rf_transcript_score", "rsem_eval_transcript_score", "transrate_transcript_score", "cats_rf_coverage_component", "cats_rf_accuracy_component", "cats_rf_local_fidelity_component", "cats_rf_integrity_component", "transrate_coverage_score", "transrate_nucleotide_score", "transrate_order_score", "transrate_segmentation_score")
metric_f_score_cor_results <- vector("list", length(metrics))
for (i in 1 : length(metrics)) {
 metric_to_analyse <- metrics[i]
 metric_f_score_cor_results[[i]] <- controlled_simulated_transcript_scores[, .("cor_coeff" = cor(get(metric_to_analyse), transcript_f_score, method = "spearman", use = "pairwise.complete.obs")), by = "assembly"]
 metric_f_score_cor_results[[i]][, "metric" := metric_to_analyse]
}
metric_f_score_cor_results <- rbindlist(metric_f_score_cor_results)

metric_f_score_cor_results[grepl("cats_rf", metric, fixed = T), "tool" := "CATS-rf"]
metric_f_score_cor_results[grepl("rsem_eval", metric, fixed = T), "tool" := "RSEM-EVAL"]
metric_f_score_cor_results[grepl("transrate", metric, fixed = T), "tool" := "TransRate"]

metric_f_score_cor_results[grepl("transcript", metric, fixed = T), "category" := "Transcript"]
metric_f_score_cor_results[grepl("coverage", metric, fixed = T), "category" := "Coverage"]
metric_f_score_cor_results[grepl("accuracy|nucleotide", metric), "category" := "Accuracy"]
metric_f_score_cor_results[grepl("local_fidelity", metric, fixed = T), "category" := "Loc. fid."]
metric_f_score_cor_results[grepl("integrity", metric, fixed = T), "category" := "Integrity"]
metric_f_score_cor_results[grepl("order", metric, fixed = T), "category" := "Order"]
metric_f_score_cor_results[grepl("segmentation", metric, fixed = T), "category" := "Segmentation"]
metric_f_score_cor_results[, "category" := factor(category, levels = c("Transcript", "Coverage", "Accuracy", "Loc. fid.", "Integrity", "Order", "Segmentation"))]

#Plotting Figure 2d (correlation between CATS-rf, RSEM-EVAL, and TransRate transcript scores and score components with transcript F-scores) (controlled simulated assemblies)
fig_2d <- ggplot(data = metric_f_score_cor_results, aes(x = category, y = cor_coeff, fill = tool)) +
 geom_boxplot(width = 0.8, size = 0.2, alpha = 0.8, color ="grey40", outlier.shape = NA) +
 theme_minimal() +
 scale_fill_manual(values = c("#E64B35", "#00A087", "#4DBBD5")) +
 theme(legend.position = "top", legend.justification = "center", legend.margin = margin(0, 0, 0, 0), legend.box.margin = margin(-3, -6, -6, -6)) +
 theme(legend.title = element_blank()) +
 theme(legend.text = element_text(size = 4)) +
 theme(legend.key.size = unit(0.22, 'cm')) +
 theme(panel.grid = element_blank()) +
 theme(panel.background = element_rect(fill = NA, color = "grey70")) +
 xlab("Score (component)") +
 ylab("Correlation with F-score") +
 scale_y_continuous(breaks = c(-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1), limits = c(-1, 1)) +
 theme(axis.title = element_text(size = 5)) +
 theme(axis.text.x = element_text(angle = 40, hjust = 0.9, size = 4.6)) +
 theme(axis.text.y = element_text(size = 4)) 

ggsave(fig_2d, file = "Figure_2d.svg", height = 1.8, width = 2.5)

#Calculating pairwise correlation between CATS-rf score components and correlations with transcript F-scores by species, coverage, mismatch rate, and assembler
get_score_f_score_corr <- function(assembly_name, subset_name) {
 if(controlled_simulated_transcript_scores_subset[, .N] != 0) {
  cor_coeffs <- sapply(metrics, function(metric) cor(controlled_simulated_transcript_scores_subset[[metric]], controlled_simulated_transcript_scores_subset[, transcript_f_score], use = "pairwise.complete.obs", method = "spearman"))
  data.table("assembly" = assembly_name, "subset" = subset_name, "metric" = names(cor_coeffs), "cor_coeff" = as.numeric(cor_coeffs))   
 }
}

pairwise_cats_rf_comp_corr_matrix_list_species_1 <- list()
pairwise_cats_rf_comp_corr_matrix_list_species_2 <- list()
pairwise_cats_rf_comp_corr_matrix_list_cov_1 <- list()
pairwise_cats_rf_comp_corr_matrix_list_cov_2 <- list()
pairwise_cats_rf_comp_corr_matrix_list_mismatch_1 <- list()
pairwise_cats_rf_comp_corr_matrix_list_mismatch_2 <- list()
pairwise_cats_rf_comp_corr_matrix_list_assembler_1 <- list()
pairwise_cats_rf_comp_corr_matrix_list_assembler_2 <- list()

cats_rf_score_f_score_corr_list_species_1 <- list()
cats_rf_score_f_score_corr_list_species_2 <- list()
cats_rf_score_f_score_corr_list_cov_1 <- list()
cats_rf_score_f_score_corr_list_cov_2 <- list()
cats_rf_score_f_score_corr_list_mismatch_1 <- list()
cats_rf_score_f_score_corr_list_mismatch_2 <- list()
cats_rf_score_f_score_corr_list_assembler_1 <- list()
cats_rf_score_f_score_corr_list_assembler_2 <- list()

for (i in 1 : length(assemblies)) {
 controlled_simulated_transcript_scores_subset <- controlled_simulated_transcript_scores[species == "S. cerevisiae" & assembly == assemblies[i]]
 pairwise_cats_rf_comp_corr_matrix_list_species_1[[i]] <- cor(controlled_simulated_transcript_scores_subset[, .(cats_rf_coverage_component, cats_rf_accuracy_component, cats_rf_local_fidelity_component, cats_rf_integrity_component)], use = "pairwise.complete.obs", method = "spearman")
 cats_rf_score_f_score_corr_list_species_1[[i]] <- get_score_f_score_corr(assemblies[i], "S. cerevisiae")
  
 controlled_simulated_transcript_scores_subset <- controlled_simulated_transcript_scores[species == "H. sapiens" & assembly == assemblies[i]]
 pairwise_cats_rf_comp_corr_matrix_list_species_2[[i]] <- cor(controlled_simulated_transcript_scores_subset[, .(cats_rf_coverage_component, cats_rf_accuracy_component, cats_rf_local_fidelity_component, cats_rf_integrity_component)], use = "pairwise.complete.obs", method = "spearman")
 cats_rf_score_f_score_corr_list_species_2[[i]] <- get_score_f_score_corr(assemblies[i], "H. sapiens")
  
 controlled_simulated_transcript_scores_subset <- controlled_simulated_transcript_scores[coverage == "1-4x" & assembly == assemblies[i]]
 pairwise_cats_rf_comp_corr_matrix_list_cov_1[[i]] <- cor(controlled_simulated_transcript_scores_subset[, .(cats_rf_coverage_component, cats_rf_accuracy_component, cats_rf_local_fidelity_component, cats_rf_integrity_component)], use = "pairwise.complete.obs", method = "spearman")
 cats_rf_score_f_score_corr_list_cov_1[[i]] <- get_score_f_score_corr(assemblies[i], "1-4x")
  
 controlled_simulated_transcript_scores_subset <- controlled_simulated_transcript_scores[coverage != "1-4x" & assembly == assemblies[i]]
 pairwise_cats_rf_comp_corr_matrix_list_cov_2[[i]] <- cor(controlled_simulated_transcript_scores_subset[, .(cats_rf_coverage_component, cats_rf_accuracy_component, cats_rf_local_fidelity_component, cats_rf_integrity_component)], use = "pairwise.complete.obs", method = "spearman")
 cats_rf_score_f_score_corr_list_cov_2[[i]] <- get_score_f_score_corr(assemblies[i], ">= 5-10x")
  
 controlled_simulated_transcript_scores_subset <- controlled_simulated_transcript_scores[mismatch == "ER 0.005" & assembly == assemblies[i]]
 pairwise_cats_rf_comp_corr_matrix_list_mismatch_1[[i]] <- cor(controlled_simulated_transcript_scores_subset[, .(cats_rf_coverage_component, cats_rf_accuracy_component, cats_rf_local_fidelity_component, cats_rf_integrity_component)], use = "pairwise.complete.obs", method = "spearman")
 cats_rf_score_f_score_corr_list_mismatch_1[[i]] <- get_score_f_score_corr(assemblies[i], "ER 0.005")
  
 controlled_simulated_transcript_scores_subset <- controlled_simulated_transcript_scores[mismatch == "ER 0.02" & assembly == assemblies[i]]
 pairwise_cats_rf_comp_corr_matrix_list_mismatch_2[[i]] <- cor(controlled_simulated_transcript_scores_subset[, .(cats_rf_coverage_component, cats_rf_accuracy_component, cats_rf_local_fidelity_component, cats_rf_integrity_component)], use = "pairwise.complete.obs", method = "spearman")
 cats_rf_score_f_score_corr_list_mismatch_2[[i]] <- get_score_f_score_corr(assemblies[i], "ER 0.02")
  
 controlled_simulated_transcript_scores_subset <- controlled_simulated_transcript_scores[assembler == "rnaSPAdes" & assembly == assemblies[i]]
 pairwise_cats_rf_comp_corr_matrix_list_assembler_1[[i]] <- cor(controlled_simulated_transcript_scores_subset[, .(cats_rf_coverage_component, cats_rf_accuracy_component, cats_rf_local_fidelity_component, cats_rf_integrity_component)], use = "pairwise.complete.obs", method = "spearman")
 cats_rf_score_f_score_corr_list_assembler_1[[i]] <- get_score_f_score_corr(assemblies[i], "rnaSPAdes")
  
 controlled_simulated_transcript_scores_subset <- controlled_simulated_transcript_scores[assembler == "Trinity" & assembly == assemblies[i]]
 pairwise_cats_rf_comp_corr_matrix_list_assembler_2[[i]] <- cor(controlled_simulated_transcript_scores_subset[, .(cats_rf_coverage_component, cats_rf_accuracy_component, cats_rf_local_fidelity_component, cats_rf_integrity_component)], use = "pairwise.complete.obs", method = "spearman")
 cats_rf_score_f_score_corr_list_assembler_2[[i]] <- get_score_f_score_corr(assemblies[i], "Trinity")
}
pairwise_cats_rf_comp_corr_matrixt_list <- list(`S. cerevisiae` = pairwise_cats_rf_comp_corr_matrix_list_species_1, `H. sapiens` = pairwise_cats_rf_comp_corr_matrix_list_species_2, `Coverage = 1-4x` = pairwise_cats_rf_comp_corr_matrix_list_cov_1, `Coverage >= 5-10x` = pairwise_cats_rf_comp_corr_matrix_list_cov_2, `ER 0.005` = pairwise_cats_rf_comp_corr_matrix_list_mismatch_1, `ER 0.02` = pairwise_cats_rf_comp_corr_matrix_list_mismatch_2, `rnaSPAdes` = pairwise_cats_rf_comp_corr_matrix_list_assembler_1, `Trinity` = pairwise_cats_rf_comp_corr_matrix_list_assembler_2)

#Plotting Supplementary Figure 1 (correlation matrices of CATS-rf score components by species, coverage, mismatch rate, and assembler) (controlled simulated assemblies)
cor_matrix_plots <- lapply(names(pairwise_cats_rf_comp_corr_matrixt_list), function(subset_category) {
 title_face <- if(subset_category %chin% c("S. cerevisiae", "H. sapiens")) "italic" else "plain"
 
 ggcorrplot(apply(simplify2array(pairwise_cats_rf_comp_corr_matrixt_list[[subset_category]]), c(1, 2), median, na.rm = T), method = "square", type = "full", ggtheme = ggplot2::theme_minimal(), colors = c("#E64B35", "white", "#4DBBD5"), show.legend = F, lab = T, show.diag = T,  lab_size = 2.1, digits = 2, tl.cex = 8, tl.srt = 0) +
  scale_x_discrete(labels = c(expression(S[c]), expression(S[a]), expression(S[l]), expression(S[i]))) +
  scale_y_discrete(labels = c(expression(S[c]), expression(S[a]), expression(S[l]), expression(S[i]))) +
  theme(axis.text = element_text(face = "italic", size = 13)) +
  ggtitle(subset_category) +
  theme(plot.title = element_text(face = title_face, size = 7))
})

supp_fig1 <- plot_grid(plotlist = cor_matrix_plots, nrow = 4, ncol = 2, align = "hv")
ggsave(supp_fig1, file = "Supplementary_Figure_1.tiff", width = 4.1, height = 6, dpi = 600, bg = "white")

cats_rf_score_f_score_corr_list <- list(`S. cerevisiae` = rbindlist(cats_rf_score_f_score_corr_list_species_1), `H. sapiens` = rbindlist(cats_rf_score_f_score_corr_list_species_2), `Coverage = 1-4x` = rbindlist(cats_rf_score_f_score_corr_list_cov_1), `Coverage >= 5-10x` = rbindlist(cats_rf_score_f_score_corr_list_cov_2), `ER 0.005` = rbindlist(cats_rf_score_f_score_corr_list_mismatch_1), `ER 0.02` = rbindlist(cats_rf_score_f_score_corr_list_mismatch_2), `rnaSPAdes` = rbindlist(cats_rf_score_f_score_corr_list_assembler_1), `Trinity` = rbindlist(cats_rf_score_f_score_corr_list_assembler_2))

#Plotting Supplementary Figure 2 (correlation between CATS-rf, RSEM-EVAL,and TransRate transcript scores and score components with transcript F-scores by species, coverage, mismatch rate, and assembler) (controlled simulated assemblies)
f_score_plots <- lapply(names(cats_rf_score_f_score_corr_list), function(subset_category) {
 cats_rf_score_f_score_corr_list_subset <- cats_rf_score_f_score_corr_list[[subset_category]]
 
 cats_rf_score_f_score_corr_list_subset[grepl("cats_rf", metric, fixed = T), "tool" := "CATS-rf"]
 cats_rf_score_f_score_corr_list_subset[grepl("rsem_eval", metric, fixed = T), "tool" := "RSEM-EVAL"]
 cats_rf_score_f_score_corr_list_subset[grepl("transrate", metric, fixed = T), "tool" := "TransRate"]
 cats_rf_score_f_score_corr_list_subset[grepl("transcript", metric, fixed = T), "category" := "Transcript"]
 cats_rf_score_f_score_corr_list_subset[grepl("coverage", metric, fixed = T), "category" := "Coverage"]
 cats_rf_score_f_score_corr_list_subset[grepl("accuracy|nucleotide", metric), "category" := "Accuracy"]
 cats_rf_score_f_score_corr_list_subset[grepl("local_fidelity", metric, fixed = T), "category" := "Loc. fid."]
 cats_rf_score_f_score_corr_list_subset[grepl("integrity", metric, fixed = T), "category" := "Integrity"]
 cats_rf_score_f_score_corr_list_subset[grepl("order", metric, fixed = T), "category" := "Order"]
 cats_rf_score_f_score_corr_list_subset[grepl("segmentation", metric, fixed = T), "category" := "Segmentation"]
 cats_rf_score_f_score_corr_list_subset[, "category" := factor(category, levels = c("Transcript", "Coverage", "Accuracy", "Loc. fid.", "Integrity", "Order", "Segmentation"))]
 
 title_face <- if(subset_category %chin% c("S. cerevisiae", "H. sapiens")) "italic" else "plain"
 ggplot(data = cats_rf_score_f_score_corr_list_subset, aes(x = category, y = cor_coeff, fill = tool)) +
  geom_boxplot(width = 0.8, size = 0.2, alpha = 0.8, color ="grey40", outlier.shape = NA) +
  theme_minimal() +
  scale_fill_manual(values = c("#E64B35", "#00A087", "#4DBBD5")) +
  theme(legend.position = "bottom", legend.justification = "center", legend.margin = margin(0, 0, 0, 0), legend.box.margin = margin(-20, -6, -6 -6)) +
  theme(legend.title = element_blank()) +
  theme(legend.text = element_text(size = 5.5)) +
  theme(legend.key.size = unit(0.4, 'cm')) +
  theme(panel.grid = element_blank()) +
  theme(panel.background = element_rect(fill = NA, color = "grey70")) +
  xlab("Score (component)") +
  ylab("Correlation with F-score") +
  scale_y_continuous(breaks = c(-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1), limits = c(-1, 1)) +
  theme(axis.title = element_text(size = 5.4)) +
  theme(axis.text.x = element_text(angle = 40, hjust = 0.9, size = 4.6)) +
  theme(axis.text.y = element_text(size = 4.3)) +
  ggtitle(subset_category) +
  theme(plot.title = element_text(size = 6.8, face = title_face))
})

supp_fig2_legend <- get_legend(f_score_plots[[1]])
supp_fig2_wo_legend <- lapply(f_score_plots, function(plot) plot + theme(legend.position = "none"))
supp_fig2_panels <- plot_grid(plotlist = supp_fig2_wo_legend, nrow = 4, ncol = 2, align = "hv")
supp_fig2 <- plot_grid(supp_fig2_panels, supp_fig2_legend, ncol = 1, rel_heights = c(1, 0.12))

ggsave(supp_fig2, file = "Supplementary_Figure_2.tiff", height = 7, width = 5, dpi = 600, bg = "white")

#Calculating correlations of CATS-rf, RSEM-EVAL, and TransRate transcript scores with transcript F-scores per species, coverage, and assembler
controlled_simulated_transcript_scores[, "sp_cov_assem" := paste(species, coverage, assembler, sep = "_")]
cor_with_f_score_assembly <- controlled_simulated_transcript_scores[, .("cats_rf_cor" = cor(cats_rf_transcript_score, transcript_f_score, use = "pairwise.complete.obs", method = "spearman"), "rsem_eval_cor" = cor(rsem_eval_transcript_score, transcript_f_score,  use = "pairwise.complete.obs", method = "spearman"), "transrate_cor" = cor(transrate_transcript_score, transcript_f_score, use = "pairwise.complete.obs", method = "spearman")), by = c("sp_cov_assem", "mismatch")]

#Creating a table with correlation results
cor_with_f_score_assembly_median_iqr <- cor_with_f_score_assembly[, .("cats_rf_cor_median" = median(cats_rf_cor), "cats_rf_cor_iqr" = IQR(cats_rf_cor), "rsem_eval_cor_median" = median(rsem_eval_cor), "rsem_eval_cor_iqr" = IQR(rsem_eval_cor), "transrate_cor_median" = median(transrate_cor), "transrate_cor_iqr" = IQR(transrate_cor)), by = "sp_cov_assem"]
cor_with_f_score_assembly_median_iqr_long <- melt(cor_with_f_score_assembly_median_iqr, id.vars = "sp_cov_assem", measure = list(cor_median = c("cats_rf_cor_median", "rsem_eval_cor_median", "transrate_cor_median"), cor_iqr = c("cats_rf_cor_iqr", "rsem_eval_cor_iqr", "transrate_cor_iqr")), variable.name = "tool")
cor_with_f_score_assembly_median_iqr_long[, "tool" := factor(tool, levels = 1 : 3, labels = c("CATS_rf", "RSEM-EVAL", "TransRate"))]

#Extracting species names, coverage, and assembler information
cor_with_f_score_assembly_median_iqr_long[, "species" := sub("_.*", "", sp_cov_assem)]
cor_with_f_score_assembly_median_iqr_long[, "coverage" := sub("^[^_]*_([^_]*)_.*", "\\1", sp_cov_assem)]
cor_with_f_score_assembly_median_iqr_long[, "assembler" := sub("^[^_]*_[^_]*_(.*)", "\\1", sp_cov_assem)]

#Adjusting factor levels for species, coverage, assembler, and tool
cor_with_f_score_assembly_median_iqr_long[, "species" := factor(species, levels = c("S. cerevisiae", "C. elegans", "D. melanogaster", "A. thaliana", "M. musculus", "H. sapiens"))]
cor_with_f_score_assembly_median_iqr_long[, "coverage" := factor(coverage, levels = rev(c("1-4x", "5-10x", "11-20x", "21-30x", "31-40x", "41-50x", "51-100x")))]
cor_with_f_score_assembly_median_iqr_long[, "assembler" := factor(assembler, levels = c("rnaSPAdes", "Trinity","IDBA-tran","SOAPdenovo-Trans"))]

#Plotting Figure 2e (correlation between CATS-rf, RSEM-EVAL, and TransRate transcript scores with transcript F-scores per species, coverage, and assembler) (controlled simulated assemblies)
fig_2e <- ggplot(data = cor_with_f_score_assembly_median_iqr_long,aes(x = coverage, y = cor_median, color = tool)) +
 geom_point(position = position_dodge(width = 0), shape = 16, size = 0.5, alpha = 0.7) +
 geom_errorbar(aes(ymin = cor_median - cor_iqr, ymax = cor_median + cor_iqr), position = position_dodge(width = 0), width = 0.2, size  = 0.2) +
 theme_minimal() +
 scale_color_manual(values = c("#E64B35", "#00A087", "#4DBBD5"), labels = c("CATS-rf", "RSEM-EVAL", "TransRate")) +
 theme(legend.position = "bottom", legend.justification = "center", legend.margin = margin(0, 0, 0, 0), legend.box.margin = margin(-6, -6, -3, -6)) +
 guides(colour = guide_legend(nrow = 1)) +
 theme(legend.title = element_blank()) +
 theme(legend.text = element_text(size = 4.5)) +
 theme(legend.key.size = unit(0.3, 'cm')) +
 facet_grid(assembler ~ species) +
 theme(strip.text.x = element_text(size = 4, face = "italic")) +
 theme(strip.text.y = element_text(size = 3)) +
 theme(panel.grid = element_blank()) +
 theme(panel.background = element_rect(fill = NA, color = "grey70")) +
 coord_flip(expand = T)  +
 xlab("Coverage") +
 ylab("Correlation with F-score") +
 scale_y_continuous(breaks = c(-1, -0.7, -0.3, 0, 0.3, 0.7, 1), labels = c(-1, -0.7, -0.3, 0, 0.3, 0.7, 1), limits = c(-1, 1)) +
 theme(axis.title = element_text(size = 4.5)) +
 theme(axis.text.x = element_text(size = 2.9)) +
 theme(axis.text.y = element_text(size = 3.3)) +
 geom_hline(aes(yintercept = 0), size = 0.25, linetype = "dashed", color = "grey70") 

ggsave(fig_2e, file = "Figure_2e.svg", height = 2.5, width = 4)

controlled_simulated_transcript_score_means <- controlled_simulated_transcript_scores[, .("cats_rf_tr_score_mean" = mean(cats_rf_transcript_score), "rsem_eval_tr_score_mean" = mean(rsem_eval_transcript_score), "transrate_tr_score_mean" = mean(transrate_transcript_score), "transcript_f_score_mean" = mean(transcript_f_score)), by = c("assembly", "species")]

#Plotting Figure 2f (correlation between CATS-rf, RSEM-EVAL, and TransRate mean transcript scores with mean transcript F-scores) (controlled simulated assemblies)
cont_sim_assem_cats_rf_f_score_cor_coeff <- sprintf("%.2f", cor(controlled_simulated_transcript_score_means[, cats_rf_tr_score_mean], controlled_simulated_transcript_score_means[, transcript_f_score_mean], method = "spearman"))
cont_sim_assem_cats_rf_f_score_dotplot <- ggplot(data = controlled_simulated_transcript_score_means, aes(x = transcript_f_score_mean, y = cats_rf_tr_score_mean, color = species)) +
 geom_point(size = 0.015, alpha = 0.8) +
 theme_minimal() +
 scale_color_npg() +
 theme(legend.position = "none") +
 theme(panel.grid = element_blank()) +
 theme(panel.background = element_rect(fill = NA, color = "grey70")) +
 xlab("Mean (F-score)") +
 ylab("Mean transcript score") +
 scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0, 1)) +
 scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0, 1)) +
 theme(axis.title = element_text(size = 5.5)) +
 theme(axis.text = element_text(size = 4.7)) +
 annotate("text", x = 0, y = 1, label = paste("r =", cont_sim_assem_cats_rf_f_score_cor_coeff), size = 1.35, hjust = 0) +
 ggtitle("CATS-rf") +
 theme(plot.title = element_text(size = 5)) 

cont_sim_assem_rsem_eval_f_score_cor_coeff <- sprintf("%.2f", cor(controlled_simulated_transcript_score_means[, rsem_eval_tr_score_mean], controlled_simulated_transcript_score_means[, transcript_f_score_mean], method = "spearman"))
cont_sim_assem_rsem_eval_f_score_dotplot <- ggplot(data = controlled_simulated_transcript_score_means, aes(x = transcript_f_score_mean, y = rsem_eval_tr_score_mean, color = species)) +
 geom_point(size = 0.015, alpha = 0.8) +
 theme_minimal() +
 scale_color_npg() +
 theme(legend.position = "none") +
 theme(panel.grid = element_blank()) +
 theme(panel.background = element_rect(fill = NA, color = "grey70")) +
 xlab("Mean (F-score)") +
 ylab("Mean transcript score") +
 scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0, 1)) +
 ylim(c(controlled_simulated_transcript_score_means[, min(rsem_eval_tr_score_mean)], controlled_simulated_transcript_score_means[, max(rsem_eval_tr_score_mean)])) +
 theme(axis.title = element_text(size = 5.5)) +
 theme(axis.text = element_text(size = 4.7)) +
 annotate("text", x = 0, y = controlled_simulated_transcript_score_means[, max(rsem_eval_tr_score_mean)], label = paste("r =", cont_sim_assem_rsem_eval_f_score_cor_coeff), size = 1.35, hjust = 0) +
 ggtitle("RSEM-EVAL") +
 theme(plot.title = element_text(size = 5)) 

cont_sim_assem_transrate_f_score_cor_coeff <- sprintf("%.2f", cor(controlled_simulated_transcript_score_means[, transrate_tr_score_mean], controlled_simulated_transcript_score_means[, transcript_f_score_mean], method = "spearman"))
cont_sim_assem_transrate_f_score_dotplot <- ggplot(data = controlled_simulated_transcript_score_means, aes(x = transcript_f_score_mean, y = transrate_tr_score_mean, color = species)) +
 geom_point(size = 0.015, alpha = 0.8) +
 theme_minimal() +
 scale_color_npg() +
 theme(legend.position = "none") +
 theme(panel.grid = element_blank()) +
 theme(panel.background = element_rect(fill = NA, color = "grey70")) +
 xlab("Mean (F-score)") +
 ylab("Mean transcript score") +
 scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0, 1)) +
 scale_y_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1), limits = c(0, 1)) +
 theme(axis.text = element_text(size = 4.7)) +
 theme(axis.title = element_text(size = 5.5)) +
 annotate("text", x = 0, y = 1, label = paste("r =", cont_sim_assem_transrate_f_score_cor_coeff), size = 1.35, hjust = 0) +
 ggtitle("TransRate") +
 theme(plot.title = element_text(size = 5)) 

fig_2f <- plot_grid(cont_sim_assem_cats_rf_f_score_dotplot, cont_sim_assem_rsem_eval_f_score_dotplot, cont_sim_assem_transrate_f_score_dotplot, nrow = 1)
ggsave(fig_2f, file = "Figure_2f.svg", height = 1.5, width = 4.6)

#Importing merged realistic simulations from file
realistically_simulated_transcript_scores <- fread("merged_realistically_simulated_transcript_scores_for_figure2.tsv")

#Extracting species names, coverage, mismatch rate, and assembler information
realistically_simulated_transcript_scores[, "assembly" := sub("realistic_", "", assembly, fixed = T)]

realistically_simulated_transcript_scores[, "species" := sub("^sim_([^_]+_[^_]+)_.*", "\\1",  assembly)]
realistically_simulated_transcript_scores[, "species" := sub("_", ". ", species, fixed = T)]
realistically_simulated_transcript_scores[, "species" := paste0(toupper(substr(species, 1, 1)), substr(species, 2, nchar(species)))]

realistically_simulated_transcript_scores[, "coverage" := sub(".*100_([^_]+)_.*", "\\1", assembly)]
realistically_simulated_transcript_scores[, "coverage" := paste0(coverage, "x")]

realistically_simulated_transcript_scores[, "assembler" := sub(".*_", "", assembly)]
realistically_simulated_transcript_scores[assembler == "RSP", "assembler" := "rnaSPAdes"]
realistically_simulated_transcript_scores[assembler == "TRI", "assembler" := "Trinity"]
realistically_simulated_transcript_scores[assembler == "IDB", "assembler" := "IDBA-tran"]
realistically_simulated_transcript_scores[assembler == "SOA", "assembler" := "SOAPdenovo-Trans"]

#Adjusting factor levels for species, coverage, mismatch, and assemblers
realistically_simulated_transcript_scores[, "species" := factor(species, levels = c("S. cerevisiae", "C. elegans", "D. melanogaster", "A. thaliana", "M. musculus", "H. sapiens"))]
realistically_simulated_transcript_scores[, "coverage" := factor(coverage, levels = c("20x", "40x", "60x", "80x"))]
realistically_simulated_transcript_scores[, "assembler" := factor(assembler, levels = c("rnaSPAdes", "Trinity", "IDBA-tran", "SOAPdenovo-Trans"))]

#Plotting Figure 2g (distribution of CATS-rf scores according to species, assembler, and coverage (realistically simulated assemblies)
fig_2g <- ggplot(data = realistically_simulated_transcript_scores, aes(x = cats_rf_transcript_score, group = coverage, color = coverage)) +
 geom_density(aes(y = ..density.. / max(..density..)), adjust = 3, linewidth = 0.25) +
 theme_minimal() +
 scale_color_manual(values = c("#4DBBD5", "#00A087", "#FDAE61", "#E64B35"), name = "Coverage") +
 theme(legend.position = "bottom", legend.justification = "center", legend.margin = margin(0, 0, 0, 0), legend.box.margin = margin(-6, -6, -3, -6)) +
 theme(legend.title = element_text(size = 5.2)) +
 theme(legend.text = element_text(size = 4)) +
 theme(legend.key.size = unit(0.25, 'cm')) + 
 facet_grid(species ~ assembler) +
 theme(strip.text.y = element_text(size = 5, face = "italic"))  +
 theme(strip.text.x = element_text(size = 5))  +
 theme(panel.grid = element_blank())  +
 theme(panel.background = element_rect(fill = NA, color = "grey70")) +
 xlab(expression(S[t])) +
 ylab("Relative density") +
 theme(axis.title = element_text(size = 6)) +
 theme(axis.title.x = element_text(face = "italic")) +
 theme(axis.text = element_text(size = 4)) 

ggsave(fig_2g, file = "Figure_2g.svg", height = 4.4, width = 4.13)

realistically_simulated_transcript_scores[, "sp_cov_assem" := paste(species, coverage, assembler, sep = "_")]
realistically_simulated_transcript_scores[, "seed" := sub(".*_([^_]*)_.*$", "\\1", assembly)]
cor_with_f_score_assembly <- realistically_simulated_transcript_scores[, .("cats_rf_cor"  = cor(cats_rf_transcript_score, transcript_f_score, use = "pairwise.complete.obs", method = "spearman"), "rsem_eval_cor" = cor(rsem_eval_transcript_score, transcript_f_score,  use = "pairwise.complete.obs", method = "spearman"), "transrate_cor" = cor(transrate_transcript_score, transcript_f_score, use = "pairwise.complete.obs", method = "spearman")), by = c("sp_cov_assem", "seed")]

cor_with_f_score_assembly_median_iqr <- cor_with_f_score_assembly[, .("cats_rf_cor_median" = median(cats_rf_cor), "cats_rf_cor_iqr" = IQR(cats_rf_cor), "rsem_eval_cor_median" = median(rsem_eval_cor), "rsem_eval_cor_iqr" = IQR(rsem_eval_cor), "transrate_cor_median" = median(transrate_cor), "transrate_cor_iqr" = IQR(transrate_cor)), by = "sp_cov_assem"]
cor_with_f_score_assembly_median_iqr_long <- melt(cor_with_f_score_assembly_median_iqr, id.vars = "sp_cov_assem", measure = list(cor_median = c("cats_rf_cor_median", "rsem_eval_cor_median", "transrate_cor_median"), cor_iqr = c("cats_rf_cor_iqr", "rsem_eval_cor_iqr", "transrate_cor_iqr")), variable.name = "tool")
cor_with_f_score_assembly_median_iqr_long[, "tool" := factor(tool, levels = 1 : 3, labels = c("CATS_rf", "RSEM-EVAL", "TransRate"))]

#Extracting species names, coverage, and assembler information
cor_with_f_score_assembly_median_iqr_long[, "species" := sub("_.*", "", sp_cov_assem)]
cor_with_f_score_assembly_median_iqr_long[, "coverage" := sub("^[^_]*_([^_]*)_.*", "\\1", sp_cov_assem)]
cor_with_f_score_assembly_median_iqr_long[, "assembler" := sub("^[^_]*_[^_]*_(.*)", "\\1", sp_cov_assem)]

#Adjusting factor levels for species, coverage, assembler, and tool
cor_with_f_score_assembly_median_iqr_long[, "species" := factor(species, levels = c("S. cerevisiae", "C. elegans", "D. melanogaster", "A. thaliana", "M. musculus", "H. sapiens"))]
cor_with_f_score_assembly_median_iqr_long[, "coverage" := factor(coverage, levels = rev(c("20x", "40x", "60x", "80x")))]
cor_with_f_score_assembly_median_iqr_long[, "assembler" := factor(assembler, levels = c("rnaSPAdes", "Trinity","IDBA-tran","SOAPdenovo-Trans"))]

#Plotting Figure 2h (correlation between CATS-rf, RSEM-EVAL, and TransRate transcript scores with transcript F-scores per species, coverage, and assembler) (realistically simulated assemblies)
fig_2h <- ggplot(data = cor_with_f_score_assembly_median_iqr_long, aes(x = coverage, y = cor_median, color = tool)) +
 geom_point(position = position_dodge(width = 0), shape = 16, size = 0.5, alpha = 0.7) +
 theme_minimal() +
 scale_color_manual(values = c("#E64B35", "#00A087", "#4DBBD5"), labels = c("CATS-rf", "RSEM-EVAL", "TransRate")) +
 theme(legend.position = "bottom", legend.justification = "center", legend.margin = margin(0, 0, 0, 0), legend.box.margin = margin(-6, -6, -3, -6)) +
 guides(colour = guide_legend(nrow = 1)) +
 theme(legend.title = element_blank()) +
 theme(legend.text = element_text(size = 4.5)) +
 theme(legend.key.size = unit(0.3, 'cm')) +
 facet_grid(assembler ~ species) +
 theme(strip.text.x = element_text(size = 4, face = "italic")) +
 theme(strip.text.y = element_text(size = 3)) +
 theme(panel.grid = element_blank()) +
 theme(panel.background = element_rect(fill = NA, color = "grey70")) +
 coord_flip(expand = T)  +
 xlab("Coverage") +
 ylab("Correlation with F-score") +
 scale_y_continuous(breaks = c(-1, -0.7, -0.3, 0, 0.3, 0.7, 1), labels = c(-1, -0.7, -0.3, 0, 0.3, 0.7, 1), limits = c(-1, 1)) +
 theme(axis.title = element_text(size = 4.5)) +
 theme(axis.text.x = element_text(size = 2.9)) +
 theme(axis.text.y = element_text(size = 3.3)) +
 geom_hline(aes(yintercept = 0), size = 0.25, linetype = "dashed", color = "grey70") 

ggsave(fig_2h, file = "fig_2h.svg", height = 2.5, width = 3.5)

realistically_simulated_transcript_score_means <- realistically_simulated_transcript_scores[, .("cats_rf_tr_score_mean" = mean(cats_rf_transcript_score), "rsem_eval_tr_score_mean" = mean(rsem_eval_transcript_score), "transrate_tr_score_mean" = mean(transrate_transcript_score), "transcript_f_score_mean" = mean(transcript_f_score)), by = c("assembly", "species")]

#Plotting Figure 2i (correlation between CATS-rf, RSEM-EVAL, and TransRate mean transcript scores with mean transcript F-scores) (realistically simulated assemblies)
real_sim_assem_cats_rf_f_score_cor_coeff <- sprintf("%.2f", cor(realistically_simulated_transcript_score_means[, cats_rf_tr_score_mean], realistically_simulated_transcript_score_means[, transcript_f_score_mean], method = "spearman"))
real_sim_assem_cats_rf_f_score_dotplot <- ggplot(data = realistically_simulated_transcript_score_means, aes(x = transcript_f_score_mean, y = cats_rf_tr_score_mean, color = species)) +
 geom_point(size = 0.015, alpha = 0.8) +
 theme_minimal() +
 scale_color_npg() +
 theme(legend.position = "none") +
 theme(panel.grid = element_blank()) +
 theme(panel.background = element_rect(fill = NA, color = "grey70")) +
 xlab("Mean (F-score)") +
 ylab("Mean transcript score") +
 scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0, 1)) +
 scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0, 1)) +
 theme(axis.title = element_text(size = 5.5)) +
 theme(axis.text = element_text(size = 4.7)) +
 annotate("text", x = 0, y = 1, label = paste("r =", real_sim_assem_cats_rf_f_score_cor_coeff), size = 1.35, hjust = 0) +
 ggtitle("CATS-rf") +
 theme(plot.title = element_text(size = 5)) 

real_sim_assem_rsem_eval_f_score_cor_coeff <- sprintf("%.2f", cor(realistically_simulated_transcript_score_means[, rsem_eval_tr_score_mean], realistically_simulated_transcript_score_means[, transcript_f_score_mean], method = "spearman"))
real_sim_assem_rsem_eval_f_score_dotplot <- ggplot(data = realistically_simulated_transcript_score_means, aes(x = transcript_f_score_mean, y = rsem_eval_tr_score_mean, color = species)) +
 geom_point(size = 0.015, alpha = 0.8) +
 theme_minimal() +
 scale_color_npg() +
 theme(legend.position = "none") +
 theme(panel.grid = element_blank()) +
 theme(panel.background = element_rect(fill = NA, color = "grey70")) +
 xlab("Mean (F-score)") +
 ylab("Mean transcript score") +
 scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0, 1)) +
 ylim(c(realistically_simulated_transcript_score_means[, min(rsem_eval_tr_score_mean)], realistically_simulated_transcript_score_means[, max(rsem_eval_tr_score_mean)])) +
 theme(axis.title = element_text(size = 5.5)) +
 theme(axis.text = element_text(size = 4.7)) +
 annotate("text", x = 0, y = realistically_simulated_transcript_score_means[, max(rsem_eval_tr_score_mean)], label = paste("r =", real_sim_assem_rsem_eval_f_score_cor_coeff), size = 1.35, hjust = 0) +
 ggtitle("RSEM-EVAL") +
 theme(plot.title = element_text(size = 5)) 

real_sim_assem_transrate_f_score_cor_coeff <- sprintf("%.2f", cor(realistically_simulated_transcript_score_means[, transrate_tr_score_mean], realistically_simulated_transcript_score_means[, transcript_f_score_mean], method = "spearman"))
real_sim_assem_transrate_f_score_dotplot <- ggplot(data = realistically_simulated_transcript_score_means, aes(x = transcript_f_score_mean, y = transrate_tr_score_mean, color = species)) +
 geom_point(size = 0.005, alpha = 0.8) +
 theme_minimal() +
 scale_color_npg() +
 theme(legend.position = "none") +
 theme(panel.grid = element_blank()) +
 theme(panel.background = element_rect(fill = NA, color = "grey70")) +
 xlab("Mean (F-score)") +
 ylab("Mean transcript score") +
 scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0, 1)) +
 scale_y_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1), limits = c(0, 1)) +
 theme(axis.text = element_text(size = 4.7)) +
 theme(axis.title = element_text(size = 5.5)) +
 annotate("text", x = 0, y = 1, label = paste("r =", real_sim_assem_transrate_f_score_cor_coeff), size = 1.35, hjust = 0) +
 ggtitle("TransRate") +
 theme(plot.title = element_text(size = 5)) 

fig_2i <- plot_grid(real_sim_assem_cats_rf_f_score_dotplot, real_sim_assem_rsem_eval_f_score_dotplot, real_sim_assem_transrate_f_score_dotplot, nrow = 1)
ggsave(fig_2i, file = "Figure_2i.svg", height = 1.5, width = 4.8)

#Importing merged public transcript scores from file
public_transcript_scores <- fread("merged_public_transcript_scores_for_figure2.tsv")

#Extracting species names, library SRA ID, and assembler
public_transcript_scores[, "species" := sub("^([^_]+_[^_]+).*", "\\1",  assembly)]
public_transcript_scores[, "species" := sub("_", ". ", species, fixed = T)]
public_transcript_scores[, "species" := paste0(toupper(substr(species, 1, 1)), substr(species, 2, nchar(species)))]

public_transcript_scores[, "library" := sub(".*(SRR[^_]*).*", "\\1", assembly)]

public_transcript_scores[, "assembler" := sub(".*_", "", assembly)]
public_transcript_scores[assembler == "RSP", "assembler" := "rnaSPAdes"]
public_transcript_scores[assembler == "TRI", "assembler" := "Trinity"]
public_transcript_scores[assembler == "IDB", "assembler" := "IDBA-tran"]
public_transcript_scores[assembler == "SOA", "assembler" := "SOAPdenovo-Trans"]

#Adjusting factor levels for species, library, and assembler
public_transcript_scores[, "species" := factor(species, levels = c("S. cerevisiae", "C. elegans", "D. melanogaster", "A. thaliana", "M. musculus", "H. sapiens"))]
public_transcript_scores[, "library" := factor(library, levels = c("SRR32108057", "SRR30685385", "SRR27822254", "SRR26147123", "SRR26001199", "SRR24356101", "SRR10815431", "SRR32732688", "SRR32058534", "SRR31107386", "SRR12868377", "SRR17736005", "SRR5224028", "SRR6266308", "SRR31530981", "SRR29606991", "SRR24876155", "SRR23933588", "SRR23870057", "SRR20326862", "SRR14560308", "SRR33339778", "SRR30855151", "SRR14160649", "SRR24948904", "SRR24636232", "SRR13159213", "SRR8422221", "SRR32903297", "SRR31588239", "SRR27369179", "SRR28908270", "SRR18855849", "SRR24134730", "SRR13981556", "SRR32139880", "SRR31785722", "SRR25646397", "SRR25732618", "SRR22548603", "SRR8357441", "SRR7741229"))]
public_transcript_scores[, "assembler" := factor(assembler, levels = c("rnaSPAdes", "Trinity", "IDBA-tran", "SOAPdenovo-Trans"))]

#Plotting Supplementary Figure 3 (distribution of CATS-rf transcript scores) (public assemblies)
supp_fig3 <- ggplot(public_transcript_scores, aes(x = cats_rf_transcript_score, color = assembler, group = assembler)) + 
 geom_density(aes(y = ..density.. / max(..density..)), adjust = 3, linewidth = 0.4) +
 theme_minimal() +
 scale_color_npg(name = "Assembler") +
 theme(legend.position = "bottom", legend.justification = "center", legend.margin=margin(0, 0, 0, 0), legend.box.margin=margin(-6, -6, -3, -6)) +
 theme(legend.title = element_text(size = 5.7)) +
 theme(legend.text = element_text(size = 5.5)) +
 theme(legend.key.size = unit(0.33, 'cm')) + 
 facet_wrap(. ~ library)  +
 theme(strip.text = element_text(size = 6.7))  +
 theme(panel.grid = element_blank())  +
 theme(panel.background = element_rect(fill = NA, color = "grey70")) +
 xlab(expression(S[t])) +
 ylab("Relative density") +
 theme(axis.title = element_text(size = 5.9)) +
 theme(axis.title.x = element_text(face = "italic")) +
 theme(axis.text = element_text(size = 5.7)) 

ggsave(supp_fig3, file = "Supplementary_Figure_3.tiff", height = 7, width = 9, dpi = 600, bg = "white")

#Calculating mean public transcript scores
public_transcript_scores <- public_transcript_scores[is.na(transcript_f_score) == F]
public_transcript_score_means <- public_transcript_scores[, .("cats_rf_transcript_score_mean" = mean(cats_rf_transcript_score), "rsem_eval_transcript_score_mean" = mean(rsem_eval_transcript_score), "transrate_transcript_score_mean" = mean(transrate_transcript_score), "transcript_f_score_mean" = mean(transcript_f_score)), by = c("assembly", "species")]

#Plotting Figure 2j (correlation between CATS-rf, RSEM-EVAL, and TransRate mean transcript scores with mean transcript F-scores) (public assemblies)
pub_tr_cats_rf_f_score_mean_cor_coeff <- sprintf("%.2f", cor(public_transcript_score_means[, cats_rf_transcript_score_mean], public_transcript_score_means[, transcript_f_score_mean], method = "spearman"))
pub_tr_cats_rf_f_score_mean_dotplot <- ggplot(data = public_transcript_score_means, aes(x = transcript_f_score_mean, y = cats_rf_transcript_score_mean, color = species)) +
 geom_point(size = 0.02, alpha = 0.8) +
 theme_minimal() +
 scale_color_npg() +
 theme(legend.position = "none") +
 theme(panel.grid = element_blank()) +
 theme(panel.background = element_rect(fill = NA, color = "grey70")) +
 xlab("Mean (F-score)") +
 ylab("Mean transcript score") +
 scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0, 1)) +
 scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0, 1)) +
 theme(axis.title = element_text(size = 5.5)) +
 theme(axis.text = element_text(size = 4.7)) +
 annotate("text", x = 0, y = 1, label = paste("r =", pub_tr_cats_rf_f_score_mean_cor_coeff), size = 1.35, hjust = 0) +
 ggtitle("CATS-rf") +
 theme(plot.title = element_text(size = 5)) 

pub_tr_rsem_eval_f_score_mean_cor_coeff <- sprintf("%.2f", cor(public_transcript_score_means[, rsem_eval_transcript_score_mean], public_transcript_score_means[, transcript_f_score_mean], method = "spearman"))
pub_tr_rsem_eval_f_score_mean_dotplot <- ggplot(data = public_transcript_score_means, aes(x = transcript_f_score_mean, y = rsem_eval_transcript_score_mean, color = species)) +
 geom_point(size = 0.02, alpha = 0.8) +
 theme_minimal() +
 scale_color_npg() +
 theme(legend.position = "none") +
 theme(panel.grid = element_blank()) +
 theme(panel.background = element_rect(fill = NA, color = "grey70")) +
 xlab("Mean (F-score)") +
 ylab("Mean transcript score") +
 scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0, 1)) +
 ylim(c(public_transcript_score_means[, min(rsem_eval_transcript_score_mean)], public_transcript_score_means[, max(rsem_eval_transcript_score_mean)])) +
 theme(axis.title = element_text(size = 5.5)) +
 theme(axis.text = element_text(size = 4.7)) +
 annotate("text", x = 0, y = public_transcript_score_means[, max(rsem_eval_transcript_score_mean)], label = paste("r =", pub_tr_rsem_eval_f_score_mean_cor_coeff), size = 1.35, hjust = 0) +
 ggtitle("RSEM-EVAL") +
 theme(plot.title = element_text(size = 5)) 

pub_tr_transrate_f_score_mean_cor_coeff <- sprintf("%.2f", cor(public_transcript_score_means[, transrate_transcript_score_mean], public_transcript_score_means[, transcript_f_score_mean], method = "spearman"))
pub_tr_transrate_f_score_mean_dotplot <- ggplot(data = public_transcript_score_means, aes(x = transcript_f_score_mean, y = transrate_transcript_score_mean, color = species)) +
 geom_point(size = 0.02, alpha = 0.8) +
 theme_minimal() +
 scale_color_npg() +
 theme(legend.position = "none") +
 theme(panel.grid = element_blank()) +
 theme(panel.background = element_rect(fill = NA, color = "grey70")) +
 xlab("Mean (F-score)") +
 ylab("Mean transcript score") +
 scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0, 1)) +
 scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0, 1)) +
 theme(axis.title = element_text(size = 5.5)) +
 theme(axis.text = element_text(size = 4.7)) +
 annotate("text", x = 0, y = 1, label = paste("r =", pub_tr_transrate_f_score_mean_cor_coeff), size = 1.35, hjust = 0) +
 ggtitle("TransRate") +
 theme(plot.title = element_text(size = 5)) 

fig_2j <- plot_grid(pub_tr_cats_rf_f_score_mean_dotplot, pub_tr_rsem_eval_f_score_mean_dotplot, pub_tr_transrate_f_score_mean_dotplot, nrow = 1)
ggsave(fig_2j, file = "Figure_2j.svg", height = 1.5, width = 4.6)