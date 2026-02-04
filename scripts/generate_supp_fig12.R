#!/usr/bin/env Rscript
#Script for generating Supplementary Figure 12
#Loading the required packages
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))

#Importing transcript-level coverage from file
transcript_coverages <- fread("coverage_table_for_supp_fig12.tsv")
setnames(transcript_coverages, c("coverage", "assembly"))
transcript_coverages[, "assembly" := sub("realistic_sim_", "", assembly, fixed = T)]
transcript_coverages[, "assembly" := sub("_coverage_table.tsv", "", assembly, fixed = T)]

#Importing the total number of reference transcripts from file
ref_tr_size_dt <- fread("ref_tr_size_for_supp_fig12.tsv")
setnames(ref_tr_size_dt, c("species", "ref_tr_N_total"))
ref_tr_size_dt[, "species" := sub("^(([^_]*_){1}[^_]*)_.*$", "\\1", species)]

#Calculating the number of unexpressed reference transcripts
lib_N_total <- transcript_coverages[, .("tr_N_total" = .N), by = "assembly"]
lib_N_total[, "species" := sub("^(([^_]*_){1}[^_]*)_.*$", "\\1", assembly)]
lib_N_total <- merge(lib_N_total, ref_tr_size_dt, by = "species")
lib_N_total[, "tr_to_add_N" := ref_tr_N_total - tr_N_total]
lib_N_total <- lib_N_total[, .(assembly, tr_to_add_N)]

#Adding completely uncovered reference transcripts
rows_to_add <- lib_N_total[, {
 new_rows <- data.table("coverage" = rep(0, tr_to_add_N), "assembly" = rep(assembly, tr_to_add_N))
 other_columns <- setdiff(colnames(transcript_coverages), colnames(new_rows))
 for (column in other_columns) {
  proto <- transcript_coverages[[column]]
  new_rows[, (column) := proto[NA_integer_]]
 }
 new_rows
}, by = "assembly"]

transcript_coverages <- rbindlist(list(transcript_coverages, rows_to_add), fill = T)[, 1 : 2]
transcript_coverages[, "seed" := as.numeric(sub(".*_", "", assembly))]
transcript_coverages[seed %in% 1:6, "replicate" := "1"]
transcript_coverages[seed %in% c(100, 200, 300, 400, 500, 600), "replicate" := "2"]
transcript_coverages[seed %in% c(1000, 2000, 3000, 4000, 5000, 6000), "replicate" := "3"]
transcript_coverages[seed %in% c(100000, 200000, 300000, 400000, 500000, 600000), "replicate" := "4"]
transcript_coverages[, "replicate" := factor(replicate)]
transcript_coverages[coverage == 0, "coverage" := 0.003]

transcript_coverages[, "species" := sub("^(([^_]*_){1}[^_]*)_.*$", "\\1", assembly)]
transcript_coverages[, "species" := sub("_", ". ", species, fixed = T)]
transcript_coverages[, "species" := paste0(toupper(substr(species, 1, 1)), substr(species, 2, nchar(species)))]
transcript_coverages[, "lib_coverage" := sub(".*100_([^_]+)_0\\.005.*", "\\1", assembly)]
transcript_coverages[, "lib_coverage" := paste0(lib_coverage, "x")]
transcript_coverages[, "species" := factor(species, levels = c("S. cerevisiae", "C. elegans", "D. melanogaster", "A. thaliana", "M. musculus", "H. sapiens"))]
transcript_coverages[, "lib_coverage" := factor(lib_coverage, levels = c("20x", "40x", "60x", "80x"))]

#Plotting Supplementary Figure 12 (transcript-level coverage distribution) (realistically simulated libraries)
supp_fig12 <- ggplot(transcript_coverages, aes(x = coverage, group = replicate, color = replicate)) +
 geom_density(aes(y = ..density.. / max(..density..)), adjust = 1.75, linewidth = 0.3) +
 theme_minimal() +
 scale_color_npg(name = "Replicate", labels = c("1", "2", "3", "4")) +
 theme(legend.position = "bottom", legend.justification = "center", legend.margin = margin(0, 0, 0, 0), legend.box.margin = margin(-6, -6, -3, -6)) +
 theme(legend.title = element_text(size = 5.2)) +
 theme(legend.text = element_text(size = 4.8)) +
 theme(legend.key.size = unit(0.3, 'cm')) + 
 facet_grid(lib_coverage ~ species) +
 theme(strip.text.x = element_text(size = 6.7, face = "italic"))  +
 theme(strip.text.y = element_text(size = 7))  +
 theme(panel.grid.major = element_line(color = "grey85", size = 0.3)) +
 theme(panel.grid.minor = element_line(color = "grey92", size = 0.2)) +
 scale_x_log10(breaks = c(0.1, 1, 10, 100, 1000, 10000), labels = c(0.1, 1, 10, 100, 1000, 10000)) +
 xlab("Transcript coverage") +
 ylab("Relative density") +
 theme(axis.title = element_text(size = 5.8)) +
 theme(axis.text.x = element_text(size = 4.2, angle = 40, hjust = 0.9)) +
 theme(axis.text.y = element_text(size = 4.2))

ggsave(supp_fig12, file = "Supplementary_Figure_12.tiff", height = 4, width = 6, dpi = 600, bg = "white")