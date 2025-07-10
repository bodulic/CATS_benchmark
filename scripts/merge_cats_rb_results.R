#!/usr/bin/env Rscript
#Script for merging CATS-rb results on simulated and public datasets (data for Figure 5)
#Loading the required package
suppressPackageStartupMessages(library(data.table))

#Importing data (function)
import_data <- function(import_regex) {
 files_in_directory <- grep(import_regex, list.files(path = "."), value = T)
 imported_dt <- list()
 for(i in 1 : length(files_in_directory)) {
  imported_dt[[i]] <- fread(files_in_directory[i], header = T)
  imported_dt[[i]][, "cats_rb_run" := files_in_directory[i]]
 }
 imported_dt <- rbindlist(imported_dt, fill = T) 
}

#Importing CATS-rb assembly scores of simulated assemblies from files
cats_rb_sim_assembly_table <- import_data(".*_CATS_rb_sim_(rel|annot)_results.*")

#Denoting metric and CATS-rb run name
cats_rb_sim_assembly_table[, "metric" := rep(c("exon_score", "transcript_score"), times = .N / 2 )]
cats_rb_sim_assembly_table[, "metric" := fifelse(grepl("_annot_", cats_rb_run, fixed = T), paste("annotation_based", metric, sep = "_"), paste("relative", metric, sep = "_"))]
cats_rb_sim_assembly_table[, "cats_rb_run" := sub("_annot_results", "", cats_rb_run)]
cats_rb_sim_assembly_table[, "cats_rb_run" := sub("_rel_results", "", cats_rb_run)]

#Setting column names
col_names <- unlist(fread("h_sapiens_CATS_rb_sim_rel_results1", header = F)[1])
colnames(cats_rb_sim_assembly_table) <- c(col_names, "cats_rb_run", "metric")

#Writing merged CATS-rb simulated assembly scores to file
write.table(cats_rb_sim_assembly_table, file = "merged_cats_rb_simulated_assembly_scores_for_figure5.tsv", sep = "\t", row.names =  F, col.names =  T, quote = F)

#Importing CATS-rb assembly scores of public assemblies from files
cats_rb_pub_assembly_table <- import_data(".*_CATS_rb_pub_(rel|annot)_results")

#Denoting metric and CATS-rb run name
cats_rb_pub_assembly_table[, "metric" := rep(c("exon_score", "transcript_score"), times = .N / 2 )]
cats_rb_pub_assembly_table[, "metric" := fifelse(grepl("_annot_", cats_rb_run, fixed = T), paste("annotation_based", metric, sep = "_"), paste("relative", metric, sep = "_"))]
cats_rb_pub_assembly_table[, "cats_rb_run" := sub("_annot_results", "", cats_rb_run)]
cats_rb_pub_assembly_table[, "cats_rb_run" := sub("_rel_results", "", cats_rb_run)]

#Writing merged CATS-rb public assembly scores to file
write.table(cats_rb_pub_assembly_table, file = "merged_cats_rb_public_assembly_scores_for_figure5.tsv", sep = "\t", row.names =  F, col.names =  T, quote = F)

#Importing CATS-rb assembly scores of simulated mutated assemblies from files
cats_rb_sim_mut_assembly_table <- import_data(".*_CATS_rb_sim_mut_rel_results")

#Denoting metric and CATS-rb run name
cats_rb_sim_mut_assembly_table[, "metric" := rep(c("relative_exon_score", "relative_transcript_score"), times = .N / 2 )]
cats_rb_sim_mut_assembly_table[, "cats_rb_run" := sub("_rel_results", "", cats_rb_run)]

#Writing merged CATS-rb simulated mutated assembly scores to file
write.table(cats_rb_sim_mut_assembly_table, file = "merged_cats_rb_simulated_mutated_assembly_scores_for_figure5.tsv", sep = "\t", row.names =  F, col.names =  T, quote = F)