#!/usr/bin/env Rscript
#Script for generating Figure 3
#Loading the required packages
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggsci))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(ggbiplot))

#Importing merged transcript scores of mutated and native assemblies from file
mut_analysis_transcript_scores <- fread("mutation_analysis_CATS_rf_transcript_scores_for_figure3.tsv")

#Extracting mutation type
mut_analysis_transcript_scores[, "assembly" := sub("controlled_", "", assembly, fixed = T)]
mut_analysis_transcript_scores[, "assembly" := sub("_1_12345", "", assembly, fixed = T)]

mut_analysis_transcript_scores[grepl("mut_ins", assembly, fixed = T), "mut_type" := "Insertion"]
mut_analysis_transcript_scores[grepl("mut_mis", assembly, fixed = T), "mut_type" := "Mismatch"]
mut_analysis_transcript_scores[grepl("mut_del", assembly, fixed = T), "mut_type" := "Deletion"]
mut_analysis_transcript_scores[grepl("mut_red", assembly, fixed = T), "mut_type" := "Redundancy"]
mut_analysis_transcript_scores[grepl("mut_frag", assembly, fixed = T), "mut_type" := "Fragmentation"]
mut_analysis_transcript_scores[grepl("mut_mult", assembly, fixed = T), "mut_type" := "Multiple"]
mut_analysis_transcript_scores[is.na(mut_type), "mut_type" := "Native"]

#Splitting assemblies with single mutation type, assemblies with multiplicative mutations, and native assemblies
mutated_transcript_scores_single <- mut_analysis_transcript_scores[mut_type %chin% c("Multiple", "Native") == F]
mutated_transcript_scores_multi <- mut_analysis_transcript_scores[mut_type == "Multiple"]
native_transcript_scores <- mut_analysis_transcript_scores[mut_type == "Native"]
rm(mut_analysis_transcript_scores)

#Extracting mutation level
mutated_transcript_scores_single[mut_type %chin% c("Insertion", "Mismatch", "Deletion") & grepl("0.15", assembly, fixed = T), "mut_level" := "1"]
mutated_transcript_scores_single[mut_type %chin% c("Insertion", "Mismatch", "Deletion") & grepl("0.3", assembly, fixed = T), "mut_level" := "2"]
mutated_transcript_scores_single[mut_type %chin% c("Insertion", "Mismatch", "Deletion") & grepl("0.45", assembly, fixed = T), "mut_level" := "3"]
mutated_transcript_scores_single[mut_type %chin% c("Insertion", "Mismatch", "Deletion") & grepl("0.6", assembly, fixed = T), "mut_level" := "4"]

mutated_transcript_scores_single[mut_type %chin% c("Redundancy", "Fragmentation") & grepl("0.2", assembly, fixed = T), "mut_level" := "1"]
mutated_transcript_scores_single[mut_type %chin% c("Redundancy", "Fragmentation") & grepl("0.4", assembly, fixed = T), "mut_level" := "2"]
mutated_transcript_scores_single[mut_type %chin% c("Redundancy", "Fragmentation") & grepl("0.6", assembly, fixed = T), "mut_level" := "3"]
mutated_transcript_scores_single[mut_type %chin% c("Redundancy", "Fragmentation") & grepl("0.8", assembly, fixed = T), "mut_level" := "4"]

mutated_transcript_scores_multi[grepl("0.1", assembly, fixed = T), "mut_level" := "1"]
mutated_transcript_scores_multi[grepl("0.2", assembly, fixed = T), "mut_level" := "2"]
mutated_transcript_scores_multi[grepl("0.3", assembly, fixed = T), "mut_level" := "3"]
mutated_transcript_scores_multi[grepl("0.4", assembly, fixed = T), "mut_level" := "4"]
mutated_transcript_scores_multi[grepl("0.5", assembly, fixed = T), "mut_level" := "5"]
mutated_transcript_scores_multi[grepl("0.6", assembly, fixed = T), "mut_level" := "6"]

native_transcript_scores[, "mut_level" := "0"]

#Adjusting factor levels for mutation type and level 
mutated_transcript_scores_single[, "mut_type" := factor(mut_type, levels = c("Insertion", "Mismatch", "Deletion", "Redundancy", "Fragmentation"))]
mutated_transcript_scores_single[, "mut_level" := factor(mut_level, levels = c("1", "2", "3", "4"))]
mutated_transcript_scores_multi[, "mut_level" := factor(mut_level, levels = c("1", "2", "3", "4", "5", "6"))]

#Calculating CATS-rf assembly scores of single-mutation and native assemblies
mutated_assembly_scores_single <- mutated_transcript_scores_single[, .(cats_rf_assembly_score = mean(cats_rf_transcript_score)), by = c("assembly", "mut_type", "mut_level")]
native_assembly_scores <- native_transcript_scores[, .(cats_rf_assembly_score = mean(cats_rf_transcript_score)), by = c("assembly", "mut_type", "mut_level")]

#Adding native assemblies to CATS-rf assembly score table (for each mutation type)
native_assembly_scores_for_insertion <- data.table(assembly = native_assembly_scores[, assembly], mut_type = "Insertion", mut_level = "0", cats_rf_assembly_score = native_assembly_scores[, cats_rf_assembly_score])
native_assembly_scores_for_mismatch <- data.table(assembly = native_assembly_scores[, assembly], mut_type = "Mismatch", mut_level = "0", cats_rf_assembly_score = native_assembly_scores[, cats_rf_assembly_score])
native_assembly_scores_for_deletion <- data.table(assembly = native_assembly_scores[, assembly], mut_type = "Deletion", mut_level = "0", cats_rf_assembly_score = native_assembly_scores[, cats_rf_assembly_score])
native_assembly_scores_for_redundancy <- data.table(assembly = native_assembly_scores[, assembly], mut_type = "Redundancy", mut_level = "0", cats_rf_assembly_score = native_assembly_scores[, cats_rf_assembly_score])
native_assembly_scores_for_fragmentation <- data.table(assembly = native_assembly_scores[, assembly], mut_type = "Fragmentation", mut_level = "0", cats_rf_assembly_score = native_assembly_scores[, cats_rf_assembly_score])

mutated_assembly_scores_single <- rbindlist(list(mutated_assembly_scores_single, native_assembly_scores_for_insertion, native_assembly_scores_for_mismatch, native_assembly_scores_for_deletion, native_assembly_scores_for_redundancy, native_assembly_scores_for_fragmentation))

#Adjusting factor levels for mutation level
mutated_assembly_scores_single[, "mut_level" := factor(mut_level, levels = c("0", "1", "2", "3", "4"))]

#Plotting Figure 3a (distribution of CATS-rf assembly scores across mutation types and levels) (single-mutation assemblies)
fig_3a <- ggplot(mutated_assembly_scores_single, aes(x = mut_level, y = cats_rf_assembly_score, fill = mut_level, color = mut_level)) +
 geom_boxplot(width = 0.8, size = 0.3, alpha = 0.75, color ="grey40", outlier.shape =  NA) +
 geom_point(position = position_jitter(width = 0.15, height = 0), size = 0.4, alpha = 0.8) +
 theme_minimal() +
 scale_fill_manual(values = c("#3C5488", "#4DBBD5", "#00A087", "#F39B7F", "#E64B35")) +
 scale_color_manual(values = c("#3C5488", "#4DBBD5", "#00A087", "#F39B7F", "#E64B35")) +
 theme(legend.position = "none")  +
 facet_wrap(. ~ mut_type, nrow = 1) +
 theme(strip.text = element_text(size = 6.9)) +
 theme(panel.grid = element_blank()) +
 theme(panel.background = element_rect(fill = NA, color = "grey70")) +
 xlab("Mutation level") +
 ylab("S") +
 ylim(c(0, 1)) +
 theme(axis.text = element_text(size = 6.5)) +
 theme(axis.title = element_text(size = 6.5)) +
 theme(axis.title.y = element_text(face = "italic"))

ggsave(fig_3a, file = "Figure_3a.svg", height = 1.8, width = 7)

#Denoting mutated transcripts (single-mutation assemblies)
mutated_transcript_scores_single[, "mutated_transcript" := fifelse(grepl("insertion_mut|mismatch_mut|deletion_mut|redundancy_mut|fragmentation_mut", transcript), "yes", "no")]
mutated_transcript_scores_single[, "mut_level_mutated" := fifelse(mutated_transcript == "yes", as.character(mut_level), "0")]

#Adjusting factor levels for mutation level
mutated_transcript_scores_single[, "mut_level_mutated" := factor(mut_level_mutated, levels = c("0", "1", "2", "3", "4"))]

#Plotting Figure 3b (distribution of coverage / accuracy / local fidelity component in mutated and native transcripts across mutation types and levels) (single-mutation assemblies)
mut_transcript_score_ins <- ggplot(mutated_transcript_scores_single[mut_type == "Insertion"], aes(x = mut_level_mutated, y = cats_rf_coverage_component, fill = mut_level_mutated)) +
 geom_boxplot(width = 0.8, size = 0.3, alpha = 0.75, color ="grey40", outlier.shape =  NA) +
 theme_minimal() +
 scale_fill_manual(values = c("#3C5488", "#4DBBD5", "#00A087", "#F39B7F", "#E64B35")) +
 theme(legend.position = "none") +
 theme(panel.grid = element_blank()) +
 theme(panel.background = element_rect(fill = NA, color = "grey70")) +
 xlab("Insertion level") +
 ylab(expression(S[c])) +
 ylim(c(0, 1)) +
 theme(axis.title = element_text(size = 6.5)) +
 theme(axis.title.y = element_text(face = "italic")) +
 theme(axis.text = element_text(size = 6.4)) +
 ggtitle("Insertion") +
 theme(plot.title = element_text(hjust = 0.5, size = 6.9)) 

mut_transcript_score_mis <- ggplot(mutated_transcript_scores_single[mut_type == "Mismatch"], aes(x = mut_level_mutated, y = cats_rf_accuracy_component, fill = mut_level_mutated)) +
 geom_boxplot(width = 0.8, size = 0.3, alpha = 0.75, color ="grey40", outlier.shape =  NA) +
 theme_minimal() +
 scale_fill_manual(values = c("#3C5488", "#4DBBD5", "#00A087", "#F39B7F", "#E64B35")) +
 theme(legend.position = "none") +
 theme(panel.grid = element_blank()) +
 theme(panel.background = element_rect(fill = NA, color = "grey70")) +
 xlab("Mismatch level") +
 ylab(expression(S[a])) +
 ylim(c(0, 1)) +
 theme(axis.title = element_text(size = 6.5)) +
 theme(axis.title.y = element_text(face = "italic")) +
 theme(axis.text = element_text(size = 6.4)) +
 ggtitle("Mismatch") +
 theme(plot.title = element_text(hjust = 0.5, size = 6.9)) 
 
mut_transcript_score_del <- ggplot(mutated_transcript_scores_single[mut_type == "Deletion"], aes(x = mut_level_mutated, y = cats_rf_local_fidelity_component, fill = mut_level_mutated)) +
 geom_boxplot(width = 0.8, size = 0.3, alpha = 0.75, color ="grey40", outlier.shape =  NA) +
 theme_minimal() +
 scale_fill_manual(values = c("#3C5488", "#4DBBD5", "#00A087", "#F39B7F", "#E64B35")) +
 theme(legend.position = "none") +
 theme(panel.grid = element_blank()) +
 theme(panel.background = element_rect(fill = NA, color = "grey70")) +
 xlab("Deletion level") +
 ylab(expression(S[l])) +
 ylim(c(0, 1)) +
 theme(axis.title = element_text(size = 6.5)) +
 theme(axis.title.y = element_text(face = "italic")) +
 theme(axis.text = element_text(size = 6.4)) +
 ggtitle("Deletion") +
 theme(plot.title = element_text(hjust = 0.5, size = 6.9)) 

fig_3b <- plot_grid(mut_transcript_score_ins, mut_transcript_score_mis, mut_transcript_score_del, nrow = 1)
ggsave(fig_3b, file = "Figure_3b.svg", height = 1.7, width = 5)
 
#Analysing CATS-rf coverage component in redundant assemblies - splitting coverage component into bins
red_transcript_cov_comps <- rbindlist(list(mutated_transcript_scores_single[mut_type == "Redundancy", .(cats_rf_coverage_component, mut_level)], native_transcript_scores[, .(cats_rf_coverage_component, mut_level)]))
red_transcript_cov_comps[, "cats_rf_coverage_component_cat" := cut(cats_rf_coverage_component, breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1))]
red_transcript_cov_comps[cats_rf_coverage_component == 0, "cats_rf_coverage_component_cat" := "0"]
red_transcript_cov_comp_N <- red_transcript_cov_comps[, .N, by = c("cats_rf_coverage_component_cat", "mut_level")]

#Adjusting factor levels for coverage bin and mutation level
red_transcript_cov_comp_N[, "cats_rf_coverage_component_cat" := factor(cats_rf_coverage_component_cat, levels = c("(0.8,1]", "(0.6,0.8]", "(0.4,0.6]", "(0.2,0.4]", "(0,0.2]", "0"))]
red_transcript_cov_comp_N[, "mut_level" := factor(mut_level, levels = c("0", "1", "2", "3", "4"))]
 
#Plotting Figure 3c (distribution of CATS-rf coverage component across redundancy levels) (single-mutation assemblies)
fig_3c <- ggplot(red_transcript_cov_comp_N, aes(x = mut_level, y = N, fill = cats_rf_coverage_component_cat)) +
 geom_bar(stat = "identity", position = "fill", width = 0.8, size = 0.3, color = "grey50", alpha = 0.8) +
 theme_minimal() +
 scale_fill_manual(values = c("#E64B35", "#F39B7F", "#B09C85", "#00A087", "#4DBBD5", "#3C5488"), name = expression(S[c]), labels = c("0.8-1", "0.6-0.8", "0.4-0.6", "0.2-0.4", "0-0.2", "0")) +
 theme(legend.position = "right", legend.margin = margin(0, 0, 0, 0)) +
 theme(legend.title = element_text(size = 5.7, face = "italic")) +
 theme(legend.key.size = unit(0.27, 'cm')) +
 theme(legend.text = element_text(size = 5))  +
 theme(panel.grid = element_blank()) +  
 theme(panel.background = element_rect(fill = NA, color = "grey70")) +
 xlab("Redundancy level") +
 ylab("% of transcripts") +
 scale_y_continuous(breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1), labels =  c(0, 20, 40, 60, 80, 100)) +
 theme(axis.title = element_text(size = 6.9)) +
 theme(axis.text = element_text(size = 6.5)) 
 
ggsave(fig_3c, file = "Figure_3c.svg", height = 1.6, width = 2.5)
 
#Calculating mean CATS-rf integrity component for fragmented assemblies
frag_transcripts_int_comp_mean <- mutated_transcript_scores_single[mut_type == "Fragmentation", .(cats_rf_integrity_component_mean = mean(cats_rf_integrity_component, na.rm = T)), by = c("assembly", "mut_level")]
native_assembly_int_comp_for_frag <- native_transcript_scores[, .(cats_rf_integrity_component_mean = mean(cats_rf_integrity_component, na.rm = T)), by = c("assembly", "mut_level")]
frag_transcripts_int_comp_mean <- rbindlist(list(frag_transcripts_int_comp_mean, native_assembly_int_comp_for_frag))

#Adjusting factor levels for mutation level
frag_transcripts_int_comp_mean[, "mut_level" := factor(mut_level, levels = c("0", "1", "2", "3", "4"))]
 
#Plotting Figure 3d (distribution of mean CATS-rf integrity component per assembly across fragmentation levels) (single-mutation assemblies)
fig_3d <- ggplot(frag_transcripts_int_comp_mean, aes(x = mut_level, y = cats_rf_integrity_component_mean, fill = mut_level, color = mut_level)) +
 geom_boxplot(width = 0.8, size = 0.3, alpha = 0.75, color ="grey40", outlier.shape =  NA) +
 geom_point(position = position_jitter(width = 0.15, height = 0), size = 0.25, alpha = 0.8) +
 theme_minimal() +
 scale_fill_manual(values = c("#3C5488", "#4DBBD5", "#00A087", "#F39B7F", "#E64B35")) +
 scale_color_manual(values = c("#3C5488", "#4DBBD5", "#00A087", "#F39B7F", "#E64B35")) +
 theme(legend.position = "none")  +
 theme(panel.grid = element_blank()) +
 theme(panel.background = element_rect(fill = NA, color = "grey70")) +
 xlab("Fragmentation level") +
 ylab(expression("Mean " * italic(S)[i])) +
 ylim(c(0, 1)) +
 theme(axis.title = element_text(size = 7)) +
 theme(axis.text = element_text(size = 6.8)) 
 
ggsave(fig_3d, file = "Figure_3d.svg", height = 1.8, width = 2)
 
#Extracting mean score components per assembly for PCA (single-mutation assemblies)
mutated_assembly_score_components <- mutated_transcript_scores_single[, .(cats_rf_coverage_component_mean = mean(cats_rf_coverage_component), cats_rf_accuracy_component_mean = mean(cats_rf_accuracy_component, na.rm = T), cats_rf_local_fidelity_component_mean = mean(cats_rf_local_fidelity_component, na.rm = T), cats_rf_integrity_component_mean = mean(cats_rf_integrity_component, na.rm = T)), by = c("assembly", "mut_type", "mut_level")]
mutated_assembly_score_components <- mutated_assembly_score_components[mut_level %chin% c("3", "4")]

#Performing PCA
mutated_assembly_score_components_pca <- prcomp(mutated_assembly_score_components[, .(cats_rf_coverage_component_mean, cats_rf_accuracy_component_mean, cats_rf_local_fidelity_component_mean, cats_rf_integrity_component_mean)], center = T, scale = F)
pca_groups <- factor(mutated_assembly_score_components[, mut_type], levels = c("Insertion", "Mismatch", "Deletion", "Redundancy", "Fragmentation"))

#Extracting % of explained variance
var_perc1 <- paste0(100 * round(mutated_assembly_score_components_pca$sdev[1], 3), "%")
var_perc2 <- paste0(100 * round(mutated_assembly_score_components_pca$sdev[2], 3), "%")

#Plotting Figure 3e (PCA of mean score components per assembly by mutation type (single-mutation assemblies, mutation levels 3 and 4))
fig_3e <- ggbiplot(mutated_assembly_score_components_pca, groups = pca_groups, ellipse = T, ellipse.level = 0.95, ellipse.linewidth = 0.5, ellipse.alpha = 0.1, circle = T, alpha = 0, var.axes = F) +
 labs(fill = "Region", color = "Region") +
 geom_point(aes(colour = mutated_assembly_score_components[, mut_type]), size = 0.5) + 
 theme_minimal() +
 scale_color_manual(values = c("#F39B7F", "#4DBBD5", "#00A087", "#E64B35", "#3C5488"), name = "Mutation type") +
 scale_fill_manual(values = c("#F39B7F", "#4DBBD5", "#00A087", "#E64B35", "#3C5488"), guide = "none") +
 theme(legend.position = "right", legend.margin = margin(0, 0, 0, 0), legend.spacing.x = unit(0.1, "cm")) +
 theme(legend.title = element_text(size = 6.2)) + 
 theme(legend.text = element_text(size = 5.5)) +
 theme(legend.key.size = unit(0.3, 'cm')) +
 theme(panel.grid = element_blank()) +
 theme(panel.background = element_rect(fill = NA, color = "grey70")) +
 xlab(paste("PC1", var_perc1)) +
 ylab(paste("PC2", var_perc2)) +
 theme(axis.title = element_text(size = 6.1)) +
 theme(axis.text = element_text(size = 6.42))

ggsave(fig_3e, file = "Figure_3e.svg", height = 1.7, width = 4.3)

#Calculating CATS-rf assembly scores of assemblies with multiplicative mutations
mutated_assembly_scores_multi <- mutated_transcript_scores_multi[, .(cats_rf_assembly_score = mean(cats_rf_transcript_score)), by = c("assembly", "mut_level")]
native_assembly_scores_for_multi <- native_transcript_scores[, .(cats_rf_assembly_score = mean(cats_rf_transcript_score, na.rm = T)), by = c("assembly", "mut_level")]
mutated_assembly_scores_multi <- rbindlist(list(mutated_assembly_scores_multi, native_assembly_scores_for_multi))

#Adjusting factor levels for mutation level
mutated_assembly_scores_multi[, "mut_level" := factor(mut_level, levels = c("0", "1", "2", "3", "4", "5", "6"))]

#Plotting Figure 3f (distribution of CATS-rf assembly scores across mutation levels (multiplicative-mutation assemblies)
fig_3f <- ggplot(mutated_assembly_scores_multi, aes(x = mut_level, y = cats_rf_assembly_score, fill = mut_level, color = mut_level)) +
 geom_boxplot(width = 0.8, size = 0.3, alpha = 0.75, color ="grey40", outlier.shape =  NA) +
 geom_point(position = position_jitter(width = 0.15, height = 0), size = 0.2, alpha = 0.8) +
 theme_minimal() +
 scale_fill_manual(values = c("#3C5488", "#4DBBD5", "#00A087", "#B09C85", "#FDAE61", "#F46D43", "#E64B35")) +
 scale_color_manual(values = c("#3C5488", "#4DBBD5", "#00A087", "#B09C85", "#FDAE61", "#F46D43", "#E64B35")) +
 theme(panel.grid = element_blank()) +
 theme(panel.background = element_rect(fill = NA, color = "grey70")) +
 theme(legend.position = "none")  +
 xlab("Mutation level") +
 ylab("S") +
 ylim(c(0, 1)) +
 theme(axis.title = element_text(size = 6.7)) +
 theme(axis.title.y = element_text(face = "italic")) +
 theme(axis.text = element_text(size = 6.36)) 

ggsave(fig_3f, file = "Figure_3f.svg", height = 1.7, width = 2.4)