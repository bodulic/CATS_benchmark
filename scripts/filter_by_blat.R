#!/usr/bin/env Rscript
#Script for filtering transcripts based on BLAT results
#Loading the required packages
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(Biostrings))

#Defining external arguments
ext_args <- commandArgs(trailingOnly = T)
THREAD_N <- as.numeric(ext_args[1])
BLAT_OUTPUT_PATH <- ext_args[2]
F_SCORE_THR <- as.numeric(ext_args[3])
QUERY_TRANSCRIPTOME_PATH <- ext_args[4]
OUT_PREF <- ext_args[5]

#Importing BLAT results from file
setDTthreads(THREAD_N)
blat_results <- fread(BLAT_OUTPUT_PATH, select = c(1, 10, 11, 15), header = F)
setnames(blat_results, c("matches", "qname", "qlen", "tlen"))

#Keeping only the best hit per transcript
setorder(blat_results, -matches)
blat_results <- unique(blat_results, by = "qname")

#Calculating F-score
blat_results[, "precision" := matches / qlen]
blat_results[, "recall" := matches / tlen]
blat_results[, "f_score" := (2 * precision * recall) / (precision + recall)]
sequences_to_keep <- blat_results[f_score >= F_SCORE_THR, qname]

#Importing query transcriptome assembly from file
query_transcriptome <- readDNAStringSet(QUERY_TRANSCRIPTOME_PATH, format = "fasta")
names(query_transcriptome) <- sub(" .*", "", names(query_transcriptome))

#Filtering transcripts
query_transcriptome_filtered <- query_transcriptome[names(query_transcriptome) %in% sequences_to_keep]

#Writing filtered query transcriptome assembly to file
writeXStringSet(query_transcriptome_filtered, filepath = paste("filtered", F_SCORE_THR, OUT_PREF, sep = "_"), format = "fasta", append = F, compress = F, compression_level = NA)