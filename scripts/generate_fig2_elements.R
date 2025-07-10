#!/usr/bin/env Rscript
#Script for generating Figure 2
#Loading the required packages
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggsci))
suppressPackageStartupMessages(library(ggcorrplot))
suppressPackageStartupMessages(library(cowplot))

#Importing merged simulated transcript scores from file
simulated_transcript_scores <- fread("merged_simulated_transcript_scores_for_figure2.tsv")

#Extracting species names, baseline coverage, mismatch rate, and assembler information
simulated_transcript_scores[, "species" := sub("^sim_([^_]+_[^_]+)_.*", "\\1",  assembly)]
simulated_transcript_scores[, "species" := sub("_", ". ", species)]
simulated_transcript_scores[, "species" := paste0(toupper(substr(species, 1, 1)), substr(species, 2, nchar(species)))]

simulated_transcript_scores[, "coverage" := sub("^.*\\.fa_100_([^_]+_[^_]+)_[^_]+_[^_]+$", "\\1", assembly)]
simulated_transcript_scores[, "coverage" := sub("_", "-", coverage)]
simulated_transcript_scores[, "coverage" := paste0(coverage, "x")]

simulated_transcript_scores[, "mismatch" := sub(".*_([0-9.]+)_[^_]+$", "\\1", assembly)]
simulated_transcript_scores[, "mismatch" := paste("ER", mismatch)]

simulated_transcript_scores[, "assembler" := sub(".*_", "", assembly)]
simulated_transcript_scores[assembler == "RSP", "assembler" := "rnaSPAdes"]
simulated_transcript_scores[assembler == "TRI", "assembler" := "Trinity"]
simulated_transcript_scores[assembler == "IDB", "assembler" := "IDBA-tran"]
simulated_transcript_scores[assembler == "SOA", "assembler" := "SOAPdenovo-Trans"]

#Calculating and writing mean transcript F-scores and CATS-rf assembly scores to file (for Figures 5e and Figure 5f)
simulated_assem_cats_rf_f_scores <- simulated_transcript_scores[, .(cats_rf_assembly_score = mean(cats_rf_transcript_score), transcript_f_score_mean = mean(transcript_f_score)), by = c("species", "coverage", "mismatch", "assembler")]
write.table(simulated_assem_cats_rf_f_scores, file = "simulated_assembly_cats_rf_f_scores_for_figure5.tsv", sep = "\t", row.names = F, col.names = T, quote = F)

#Adjusting factor levels for species, coverage, mismatch, and assemblers
simulated_transcript_scores[, "species" := factor(species, levels = c("S. cerevisiae", "C. elegans", "D. melanogaster", "A. thaliana", "M. musculus", "H. sapiens"))]
simulated_transcript_scores[, "coverage" := factor(coverage, levels = c("1-4x", "5-10x", "11-20x", "21-30x", "31-40x", "41-50x", "51-100x"))]
simulated_transcript_scores[, "mismatch" := factor(mismatch, levels = c("ER 0.005", "ER 0.01", "ER 0.02"))]
simulated_transcript_scores[, "assembler" := factor(assembler, levels = c("rnaSPAdes", "Trinity", "IDBA-tran", "SOAPdenovo-Trans"))]

#Plotting Figure 2a (distribution of CATS-rf scores according to species and assembler) (simulated assemblies)
fig_2a <- ggplot(data = simulated_transcript_scores, aes(x = cats_rf_transcript_score, group = assembler, color = assembler)) +
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

#Plotting Figure 2b (distribution of CATS-rf scores according to coverage, mismatch rate, and assembler) (simulated assemblies)
fig_2b <- ggplot(data = simulated_transcript_scores, aes(x = cats_rf_transcript_score, group = coverage, color = coverage)) +
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
pairwise_cats_rf_comp_corr_matrix <- cor(simulated_transcript_scores[, .(cats_rf_coverage_component, cats_rf_accuracy_component, cats_rf_local_fidelity_component, cats_rf_integrity_component)], use = "pairwise.complete.obs", method = "spearman")

#Plotting Figure 2c (correlation matrix of CATS-rf score components) (simulated assemblies)
fig_2c <- ggcorrplot(pairwise_cats_rf_comp_corr_matrix, method = "square", type = "full", ggtheme = ggplot2::theme_minimal(), colors = c("#E64B35", "white", "#4DBBD5"), show.legend = F, lab = T, show.diag = T, lab_size = 2.1, digits = 2, tl.cex = 8, tl.srt = 0) +
 scale_x_discrete(labels = c(expression(S[c]), expression(S[a]), expression(S[l]), expression(S[i]))) + 
 scale_y_discrete(labels = c(expression(S[c]), expression(S[a]), expression(S[l]), expression(S[i]))) +
 theme(axis.text = element_text(face = "italic", size = 13))

ggsave(fig_2c, file = "figure_2c.svg", height = 1.6, width = 1.6)

#Calculating correlations of CATS-rf, RSEM-EVAL, and TransRate transcript scores and score components with transcript F-scores
metrics <- c("cats_rf_transcript_score", "rsem_eval_transcript_score", "transrate_transcript_score", "cats_rf_coverage_component", "cats_rf_accuracy_component", "cats_rf_local_fidelity_component", "cats_rf_integrity_component", "transrate_coverage_score", "transrate_nucleotide_score", "transrate_order_score", "transrate_segmentation_score")
cor_results <- c()

for (i in 1 : length(metrics)) {
 metric_to_analyse <- metrics[i]
 cor_coeff <- cor(simulated_transcript_scores[, ..metric_to_analyse], simulated_transcript_scores[, transcript_f_score], use = "pairwise.complete.obs",  method = "spearman")
 cor_results[i] <- cor_coeff
}

#Creating a table with correlation results
cor_results_table <- data.table(metric = metrics, corr_coeff = cor_results)
cor_results_table[grepl("cats_rf", metric), "tool" := "CATS-rf"]
cor_results_table[grepl("rsem_eval", metric), "tool" := "RSEM-EVAL"]
cor_results_table[grepl("transrate", metric), "tool" := "TransRate"]

#Adjusting factor levels for scores / score components
cor_results_table[, "category" := c("Transcript", "Transcript", "Transcript", "Coverage", "Accuracy", "Loc. fid.", "Integrity", "Coverage", "Accuracy", "Order", "Segmentation")]
cor_results_table[, "category" := factor(category, levels = c("Transcript", "Coverage", "Accuracy", "Loc. fid.", "Integrity", "Order", "Segmentation"))]

#Plotting Figure 2d (correlation between CATS-rf, RSEM-EVAL, and TransRate transcript scores and score components with transcript F-scores) (simulated assemblies)
fig_2d <- ggplot(data = cor_results_table, aes(x = category, y = corr_coeff, fill = tool)) +
 geom_bar(stat = "identity", position = "dodge", width = 0.75, size = 0.3, color = "grey50", alpha = 0.85) +
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
 scale_y_continuous(breaks = c(-0.1, 0, 0.25, 0.5, 0.75, 1), limits = c(0, 1)) +
 theme(axis.title = element_text(size = 5)) +
 theme(axis.text.x = element_text(angle = 40, hjust = 0.9, size = 4.6)) +
 theme(axis.text.y = element_text(size = 4)) 

ggsave(fig_2d, file = "Figure_2d.svg", height = 1.8, width = 2.3)
 
#Calculating correlations of CATS-rf, RSEM-EVAL, and TransRate transcript scores with transcript F-scores per species, coverage, and assembler
simulated_transcript_scores[, "sp_cov_assem" := paste(species, coverage, assembler, sep = "_")]

all_sp_cov_assem <- simulated_transcript_scores[, unique(sp_cov_assem)]
cats_rf_cor_coeffs <- c()
rsem_eval_cor_coeffs <- c()
transrate_cor_coeffs <- c()

for (i in 1 : length(all_sp_cov_assem)) {
 simulated_transcript_scores_filtered <- simulated_transcript_scores[sp_cov_assem == all_sp_cov_assem[i]]
 cats_rf_cor_coeffs[i] <- cor(simulated_transcript_scores_filtered[, cats_rf_transcript_score], simulated_transcript_scores_filtered[, transcript_f_score], use = "pairwise.complete.obs", method = "spearman")
 rsem_eval_cor_coeffs[i] <- cor(simulated_transcript_scores_filtered[, rsem_eval_transcript_score], simulated_transcript_scores_filtered[, transcript_f_score], use = "pairwise.complete.obs", method = "spearman")
 transrate_cor_coeffs[i] <- cor(simulated_transcript_scores_filtered[, transrate_transcript_score], simulated_transcript_scores_filtered[, transcript_f_score], use = "pairwise.complete.obs", method = "spearman")
}

#Creating a table with correlation results
sp_cov_assem_cor_table <- data.table(sp_cov_assem = all_sp_cov_assem, cats_rf_cor_coeff = cats_rf_cor_coeffs, rsem_eval_cor_coeff = rsem_eval_cor_coeffs, transrate_cor_coeff = transrate_cor_coeffs)
sp_cov_assem_cor_table_longer <- melt(sp_cov_assem_cor_table, measure.vars = c("cats_rf_cor_coeff", "rsem_eval_cor_coeff", "transrate_cor_coeff"),  variable.name = "tool", value.name = "cor_coeff")

#Extracting species names, baseline coverage, and assembler information
sp_cov_assem_cor_table_longer[, "species" := sub("_.*", "", sp_cov_assem)]
sp_cov_assem_cor_table_longer[, "coverage" := sub("^[^_]*_([^_]*)_.*", "\\1", sp_cov_assem)]
sp_cov_assem_cor_table_longer[, "assembler" := sub("^[^_]*_[^_]*_(.*)", "\\1", sp_cov_assem)]

#Adjusting factor levels for species, coverage, assembler, and tool
sp_cov_assem_cor_table_longer[, "species" := factor(species, levels = c("S. cerevisiae", "C. elegans", "D. melanogaster", "A. thaliana", "M. musculus", "H. sapiens"))]
sp_cov_assem_cor_table_longer[, "coverage" := factor(coverage, levels = rev(c("1-4x", "5-10x", "11-20x", "21-30x", "31-40x", "41-50x", "51-100x")))]
sp_cov_assem_cor_table_longer[, "assembler" := factor(assembler, levels = c("rnaSPAdes", "Trinity","IDBA-tran","SOAPdenovo-Trans"))]
sp_cov_assem_cor_table_longer[, "tool" := factor(tool, levels = c("cats_rf_cor_coeff", "rsem_eval_cor_coeff", "transrate_cor_coeff"))]

#Plotting Figure 2e (correlation between CATS-rf, RSEM-EVAL, and TransRate transcript scores with transcript F-scores per species, coverage, and assembler) (simulated assemblies)
fig_2e <- ggplot(data = sp_cov_assem_cor_table_longer, aes(x = coverage, y = cor_coeff, color = tool)) +
 geom_point(shape = 16, size = 0.5, alpha = 0.9) +
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
 coord_flip()  +
 xlab("Coverage") +
 ylab("Correlation with F-score") +
 scale_y_continuous(breaks = c(-1, -0.7,-0.3, 0, 0.3, 0.7, 1), labels = c(-1, -0.7, -0.3, 0, 0.3, 0.7, 1), limits = c(-1, 1)) +
 theme(axis.title = element_text(size = 4.5)) +
 theme(axis.text.x = element_text(size = 2.9)) +
 theme(axis.text.y = element_text(size = 3.3)) +
 geom_hline(aes(yintercept = 0), size = 0.25, linetype = "dashed", color = "grey70") 

ggsave(fig_2e, file = "Figure_2e.svg", height = 2.5, width = 4)

#Importing merged simulated assembly scores from file
simulated_assembly_scores <- fread("merged_simulated_assembly_scores_for_figure2.tsv")

#Extracting species names
simulated_assembly_scores[, "species" := sub("^sim_([^_]+_[^_]+)_.*", "\\1",  assembly)]
simulated_assembly_scores[, "species" := sub("_", ". ", species)]
simulated_assembly_scores[, "species" := paste0(toupper(substr(species, 1, 1)), substr(species, 2, nchar(species)))]

#Adjusting factor levels for species
simulated_assembly_scores[, "species" := factor(species, levels = c("S. cerevisiae", "C. elegans", "D. melanogaster", "A. thaliana", "M. musculus", "H. sapiens"))]

#Plotting Figure 2f (correlation between CATS-rf, RSEM-EVAL, and TransRate assembly scores with mean transcript F-scores) (simulated assemblies)
sim_assem_cats_rf_f_score_cor_coeff <- sprintf("%.2f", cor(simulated_assembly_scores[, cats_rf_assembly_score], simulated_assembly_scores[, transcript_f_score_mean], method = "spearman"))
sim_assem_cats_rf_f_score_dotplot <- ggplot(data = simulated_assembly_scores, aes(x = transcript_f_score_mean, y = cats_rf_assembly_score, color = species)) +
 geom_point(size = 0.015, alpha = 0.8) +
 theme_minimal() +
 scale_color_npg() +
 theme(legend.position = "none") +
 theme(panel.grid = element_blank()) +
 theme(panel.background = element_rect(fill = NA, color = "grey70")) +
 xlab("Mean (F-score)") +
 ylab("Assembly score") +
 scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0, 1)) +
 scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0, 1)) +
 theme(axis.title = element_text(size = 5.5)) +
 theme(axis.text = element_text(size = 4.7)) +
 annotate("text", x = 0, y = 1, label = paste("r =", sim_assem_cats_rf_f_score_cor_coeff), size = 1.35, hjust = 0) +
 ggtitle("CATS-rf") +
 theme(plot.title = element_text(size = 5)) 

sim_assem_rsem_eval_f_score_cor_coeff <- sprintf("%.2f", cor(simulated_assembly_scores[, rsem_eval_assembly_score], simulated_assembly_scores[, transcript_f_score_mean], method = "spearman"))
sim_assem_rsem_eval_f_score_dotplot <- ggplot(data = simulated_assembly_scores, aes(x = transcript_f_score_mean, y = rsem_eval_assembly_score, color = species)) +
 geom_point(size = 0.015, alpha = 0.8) +
 theme_minimal() +
 scale_color_npg() +
 theme(legend.position = "none") +
 theme(panel.grid = element_blank()) +
 theme(panel.background = element_rect(fill = NA, color = "grey70")) +
 xlab("Mean (F-score)") +
 ylab("Assembly score") +
 scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0, 1)) +
 ylim(c(simulated_assembly_scores[, min(rsem_eval_assembly_score)], simulated_assembly_scores[, max(rsem_eval_assembly_score)])) +
 theme(axis.title = element_text(size = 5.5)) +
 theme(axis.text = element_text(size = 4.7)) +
 annotate("text", x = 0, y = simulated_assembly_scores[, max(rsem_eval_assembly_score)], label = paste("r =", sim_assem_rsem_eval_f_score_cor_coeff), size = 1.35, hjust = 0) +
 ggtitle("RSEM-EVAL") +
 theme(plot.title = element_text(size = 5)) 

sim_assem_transrate_f_score_cor_coeff <- sprintf("%.2f", cor(simulated_assembly_scores[, transrate_assembly_score], simulated_assembly_scores[, transcript_f_score_mean], method = "spearman"))
sim_assem_transrate_f_score_dotplot <- ggplot(data = simulated_assembly_scores, aes(x = transcript_f_score_mean, y = transrate_assembly_score, color = species)) +
 geom_point(size = 0.015, alpha = 0.8) +
 theme_minimal() +
 scale_color_npg() +
 theme(legend.position = "none") +
 theme(panel.grid = element_blank()) +
 theme(panel.background = element_rect(fill = NA, color = "grey70")) +
 xlab("Mean (F-score)") +
 ylab("Assembly score") +
 scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0, 1)) +
 scale_y_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1), limits = c(0, 1)) +
 theme(axis.text = element_text(size = 4.7)) +
 theme(axis.title = element_text(size = 5.5)) +
 annotate("text", x = 0, y = 1, label = paste("r =", sim_assem_transrate_f_score_cor_coeff), size = 1.35, hjust = 0) +
 ggtitle("TransRate") +
 theme(plot.title = element_text(size = 5)) 

fig_2f <- plot_grid(sim_assem_cats_rf_f_score_dotplot, sim_assem_rsem_eval_f_score_dotplot, sim_assem_transrate_f_score_dotplot, nrow = 1)
ggsave(fig_2f, file = "Figure_2f.svg", height = 1.5, width = 4.6)

#Importing merged public transcript scores from file
public_transcript_scores <- fread("merged_public_transcript_scores_for_figure2.tsv")

#Extracting species names, library SRA ID, and assembler
public_transcript_scores[, "species" := sub("^([^_]+_[^_]+).*", "\\1",  assembly)]
public_transcript_scores[, "species" := sub("_", ". ", species)]
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

#Plotting Extended data Figure 1 (distribution of CATS-rf transcript scores) (public assemblies)
fig_ext1 <- ggplot(public_transcript_scores, aes(x = cats_rf_transcript_score, color = assembler, group = assembler)) + 
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

ggsave(fig_ext1, file = "Extended_data_figure_1.tiff", height = 7, width = 9, dpi = 600, bg = "white")

#Calculating mean public transcript scores
public_transcript_scores <- public_transcript_scores[is.na(transcript_f_score) == F]
public_transcript_score_means <- public_transcript_scores[, .(cats_rf_transcript_score_mean = mean(cats_rf_transcript_score), rsem_eval_transcript_score_mean = mean(rsem_eval_transcript_score), transrate_transcript_score_mean = mean(transrate_transcript_score), transcript_f_score_mean = mean(transcript_f_score)), by = c("assembly", "species")]

#Plotting Figure 2g (correlation between CATS-rf, RSEM-EVAL, and TransRate mean transcript scores with mean transcript F-scores) (public assemblies)
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

fig_2g <- plot_grid(pub_tr_cats_rf_f_score_mean_dotplot, pub_tr_rsem_eval_f_score_mean_dotplot, pub_tr_transrate_f_score_mean_dotplot, nrow = 1)
ggsave(fig_2g, file = "Figure_2g.svg", height = 1.5, width = 4.6)