#!/usr/bin/env Rscript
#F-score calculation script
#Loading the required packages
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(Biostrings))

#Defining external arguments
ext_args <- commandArgs(trailingOnly = T)
THREAD_N <- as.numeric(ext_args[1])
CRB_OUTPUT_PATH <- ext_args[2]
QUERY_TRANSCRIPTOME_PATH <- ext_args[3]
OUT_PREF <- ext_args[4]

#Importing CRB-BLAST results from file
setDTthreads(THREAD_N)
crb_results <- fread(CRB_OUTPUT_PATH, header = F, select = c(1, 2, 3, 4, 9, 10))
setnames(crb_results, c("qname", "tname", "perc_id", "al_len", "qlen", "tlen"))

#Calculating F-score
crb_results[, "matches" := al_len * (perc_id / 100)]
crb_results[, "precision" := matches / qlen]
crb_results[, "recall" := matches / tlen]
crb_results[, "f_score" := (2 * precision * recall) / (precision + recall)]

#Importing query transcriptome assembly from file
query_transcriptome <- readDNAStringSet(QUERY_TRANSCRIPTOME_PATH, format = "fasta")
names(query_transcriptome) <- sub(" .*", "", names(query_transcriptome))

#Adding sequences without hits (F-score = 0)
missing_query_transcripts <- query_transcriptome[names(query_transcriptome) %in% crb_results[, qname] == F]
crb_results <- rbindlist(list(crb_results, data.table(qname = names(missing_query_transcripts), tname = NA, perc_id = NA, al_len = 0, qlen = width(missing_query_transcripts), tlen = NA, matches = 0, precision = 0, recall = 0, f_score = 0)))

#Writing F-scores to file
write.table(crb_results, file = paste(OUT_PREF, "f_scores", sep = "_"), sep = "\t", row.names = F, col.names = T, quote = F)