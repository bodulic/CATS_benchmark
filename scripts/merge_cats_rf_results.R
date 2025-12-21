#!/usr/bin/env Rscript
#Script for merging CATS-rf results from controlled simulated, realistically simulated, and public datasets (data for Figure 2)
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

#Importing CATS-rf transcript scores of controlled simulated assemblies from files
cats_rf_cont_sim_transcript_table <- import_data("controlled_sim_.*_CATS_rf_transcript_scores.tsv")
setnames(cats_rf_cont_sim_transcript_table, c("transcript", "cats_rf_coverage_component", "cats_rf_accuracy_component", "cats_rf_local_fidelity_component", "cats_rf_integrity_component", "cats_rf_transcript_score", "assembly"))

#Importing RSEM-EVAL transcript scores of controlled simulated assemblies from files
rsem_eval_cont_sim_transcript_table <- import_data("controlled_sim_.*_rsem_eval_transcript_scores")
setnames(rsem_eval_cont_sim_transcript_table, c("transcript", "rsem_eval_transcript_score", "assembly"))

#Importing TransRate transcript scores of controlled simulated assemblies from files
transrate_cont_sim_transcript_table <- import_data("controlled_sim_.*_transrate_transcript_scores")
setnames(transrate_cont_sim_transcript_table, c("transcript", "transrate_coverage_score", "transrate_nucleotide_score", "transrate_order_score", "transrate_segmentation_score", "transrate_transcript_score", "assembly"))

#Importing transcript F-scores of controlled simulated assemblies from files
f_score_cont_sim_transcript_table <- import_data("controlled_sim_.*_f_scores")
setnames(f_score_cont_sim_transcript_table, c("transcript", "transcript_f_score", "assembly"))

#Merging CATS-rf, RSEM-EVAL and TransRate controlled simulated transcript scores with F-scores
merged_cont_sim_transcript_table <- merge(cats_rf_cont_sim_transcript_table, rsem_eval_cont_sim_transcript_table, by = c("assembly", "transcript"))
merged_cont_sim_transcript_table <- merge(merged_cont_sim_transcript_table, transrate_cont_sim_transcript_table, by = c("assembly", "transcript"))
merged_cont_sim_transcript_table <- merge(merged_cont_sim_transcript_table, f_score_cont_sim_transcript_table, by = c("assembly", "transcript"))

#Writing merged controlled simulated transcript scores to file
write.table(merged_cont_sim_transcript_table, file = "merged_controlled_simulated_transcript_scores_for_figure2.tsv", sep = "\t", row.names = F, col.names = T, quote = F)

#Importing CATS-rf transcript scores of realistically simulated assemblies from files
cats_rf_real_sim_transcript_table <- import_data("realistic_sim_.*_CATS_rf_transcript_scores.tsv")
setnames(cats_rf_real_sim_transcript_table, c("transcript", "cats_rf_coverage_component", "cats_rf_accuracy_component", "cats_rf_local_fidelity_component", "cats_rf_integrity_component", "cats_rf_transcript_score", "assembly"))

#Importing RSEM-EVAL transcript scores of realistically simulated assemblies from files
rsem_eval_real_sim_transcript_table <- import_data("realistic_sim_.*_rsem_eval_transcript_scores")
setnames(rsem_eval_real_sim_transcript_table, c("transcript", "rsem_eval_transcript_score", "assembly"))

#Importing TransRate transcript scores of realistically simulated assemblies from files
transrate_real_sim_transcript_table <- import_data("realistic_sim_.*_transrate_transcript_scores")
setnames(transrate_real_sim_transcript_table, c("transcript", "transrate_coverage_score", "transrate_nucleotide_score", "transrate_order_score", "transrate_segmentation_score", "transrate_transcript_score", "assembly"))

#Importing transcript F-scores of realistically simulated assemblies from files
f_score_real_sim_transcript_table <- import_data("realistic_sim_.*_f_scores")
setnames(f_score_real_sim_transcript_table, c("transcript", "transcript_f_score", "assembly"))

#Merging CATS-rf, RSEM-EVAL and TransRate realistically simulated transcript scores with F-scores
merged_real_sim_transcript_table <- merge(cats_rf_real_sim_transcript_table, rsem_eval_real_sim_transcript_table, by = c("assembly", "transcript"))
merged_real_sim_transcript_table <- merge(merged_real_sim_transcript_table, transrate_real_sim_transcript_table, by = c("assembly", "transcript"))
merged_real_sim_transcript_table <- merge(merged_real_sim_transcript_table, f_score_real_sim_transcript_table, by = c("assembly", "transcript"))

#Writing merged realistically simulated transcript scores to file
write.table(merged_real_sim_transcript_table, file = "merged_realistically_simulated_transcript_scores_for_figure2.tsv", sep = "\t", row.names = F, col.names = T, quote = F)

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