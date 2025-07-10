#!/usr/bin/env Rscript
#Script for merging CATS-rf results on simulated and public datasets (data for Figure 2)
#Loading the required package
suppressPackageStartupMessages(library(data.table))

#Importing data (function)
import_data <- function(import_regex) {
 files_in_directory <- grep(import_regex, list.files(path = "."), value = T)
 imported_dt <- list()
 for(i in 1 : length(files_in_directory)) {
  imported_dt[[i]] <- fread(files_in_directory[i], header = F)
  imported_dt[[i]][, "filename" := files_in_directory[i]]
 }
 imported_dt <- rbindlist(imported_dt) 
 remove_regex <- sub(".*\\*", "", import_regex)
 imported_dt[, "filename" := sub(remove_regex, "", filename)]
}

#Importing CATS-rf transcript scores of simulated assemblies from files
cats_rf_sim_transcript_table <- import_data("sim_.*_CATS_rf_transcript_scores.tsv")
setnames(cats_rf_sim_transcript_table, c("transcript", "cats_rf_coverage_component", "cats_rf_accuracy_component", "cats_rf_local_fidelity_component", "cats_rf_integrity_component", "cats_rf_transcript_score", "assembly"))

#Importing RSEM-EVAL transcript scores of simulated assemblies from files
rsem_eval_sim_transcript_table <- import_data("sim_.*_rsem_eval_transcript_scores")
setnames(rsem_eval_sim_transcript_table, c("transcript", "rsem_eval_transcript_score", "assembly"))

#Importing TransRate transcript scores of simulated assemblies from files
transrate_sim_transcript_table <- import_data("sim_.*_transrate_transcript_scores")
setnames(transrate_sim_transcript_table, c("transcript", "transrate_coverage_score", "transrate_nucleotide_score", "transrate_order_score", "transrate_segmentation_score", "transrate_transcript_score", "assembly"))

#Importing transcript F-scores of simulated assemblies from files
f_score_sim_transcript_table <- import_data("sim_.*_f_scores")
setnames(f_score_sim_transcript_table, c("transcript", "transcript_f_score", "assembly"))

#Merging CATS-rf, RSEM-EVAL and TransRate simulated transcript scores with F-scores
merged_sim_transcript_table <- merge(cats_rf_sim_transcript_table, rsem_eval_sim_transcript_table, by = c("assembly", "transcript"))
merged_sim_transcript_table <- merge(merged_sim_transcript_table, transrate_sim_transcript_table, by = c("assembly", "transcript"))
merged_sim_transcript_table <- merge(merged_sim_transcript_table, f_score_sim_transcript_table, by = c("assembly", "transcript"))

#Writing merged simulated transcript scores to file
write.table(merged_sim_transcript_table, file = "merged_simulated_transcript_scores_for_figure2.tsv", sep = "\t", row.names = F, col.names = T, quote = F)

#Calculating CATS-rf simulated assembly scores
cats_rf_sim_assembly_table <- cats_rf_sim_transcript_table[, .(cats_rf_assembly_score = mean(transcript_score)), by = "assembly"]

#Importing RSEM-EVAL simulated assembly scores from files
rsem_eval_sim_assembly_table <- import_data("sim_.*_rsem_eval_assembly_score")
setnames(rsem_eval_sim_assembly_table, c("rsem_eval_assembly_score", "assembly"))

#Importing TransRate simulated assembly scores from files
transrate_sim_assembly_table <- import_data("sim_.*_transrate_assembly_score")
setnames(transrate_sim_assembly_table, c("transrate_assembly_score", "assembly"))

#Calculating mean F-scores of simulated transcripts
f_score_sim_transcript_mean_table <- f_score_sim_transcript_table[, .(transcript_f_score_mean = mean(transcript_f_score)), by = "assembly"]

#Merging CATS-rf, RSEM-EVAL, and TransRate simulated assembly scores with mean F-scores
merged_sim_assembly_table <- merge(cats_rf_sim_assembly_table, rsem_eval_sim_assembly_table, by = "assembly")
merged_sim_assembly_table <- merge(merged_sim_assembly_table, transrate_sim_assembly_table, by = "assembly")
merged_sim_assembly_table <- merge(merged_sim_assembly_table, f_score_sim_transcript_mean_table, by = "assembly")

#Writing merged simulated assembly scores to file
write.table(merged_sim_assembly_table, file = "merged_simulated_assembly_scores_for_figure2.tsv", sep = "\t", row.names = F, col.names = T, quote = F)

#Importing CATS-rf transcript scores of public assemblies from files
cats_rf_pub_transcript_table <- import_data("pub_.*_CATS_rf_transcript_scores.tsv")
setnames(cats_rf_pub_transcript_table, c("transcript", "cats_rf_coverage_component", "cats_rf_accuracy_component", "cats_rf_local_fidelity_component", "cats_rf_integrity_component", "cats_rf_transcript_score", "assembly"))

#Importing RSEM-EVAL transcript scores of public assemblies from files
rsem_eval_pub_transcript_table <- import_data("pub_.*_rsem_eval_transcript_scores")
setnames(rsem_eval_pub_transcript_table, c("transcript", "rsem_eval_transcript_score", "assembly"))

#Importing TransRate transcript scores of public assemblies from files
transrate_pub_transcript_table <- import_data("pub_.*_transrate_transcript_scores")
setnames(transrate_pub_transcript_table, c("transcript", "transrate_coverage_score", "transrate_nucleotide_score", "transrate_order_score", "transrate_segmentation_score", "transrate_transcript_score", "assembly"))

#Importing transcript F-scores of public assemblies from files
f_score_pub_transcript_table <- import_data("pub_.*_f_scores")
setnames(f_score_pub_transcript_table, c("transcript", "transcript_f_score", "assembly"))

#Merging CATS-rf, RSEM-EVAL and TransRate public transcript scores with F-score
merged_pub_transcript_table <- merge(cats_rf_pub_transcript_table, rsem_eval_pub_transcript_table, by = c("assembly", "transcript"))
merged_pub_transcript_table <- merge(merged_pub_transcript_table, transrate_pub_transcript_table, by = c("assembly", "transcript"))
merged_pub_transcript_table <- merge(merged_pub_transcript_table, f_score_pub_transcript_table, by = c("assembly", "transcript"), all.x = T)

#Writing merged public transcript scores to file
write.table(merged_pub_transcript_table, file = "merged_public_transcript_scores_for_figure2.tsv", sep = "\t", row.names = F, col.names = T, quote = F)