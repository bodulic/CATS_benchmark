#!/usr/bin/env Rscript
#Read simulation script
#Loading the required packages
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(polyester))
options(scipen = 999)

#Defining external arguments
ext_args <- commandArgs(trailingOnly = T)
TRANSCRIPTOME_PATH <- ext_args[1]
SIM_TYPE <- ext_args[2]
COVERAGE <- ext_args[3]
SEED <- as.numeric(ext_args[4])
EXPRESSED_PROP_RANGE <- ext_args[5]
READ_LENGTH <- as.numeric(ext_args[6])
MISMATCH_RATE <- as.numeric(ext_args[7])
OUT_PREF <- ext_args[8]

#Importing transcriptome assembly from file
transcriptome <- readDNAStringSet(TRANSCRIPTOME_PATH, format = "fasta")

#Determining the number of simulated reads
if(SIM_TYPE == "controlled") {
 expressed_tr_prop_sample <- "1"   
 coverage_lower_boundary <- as.numeric(sub("_.*", "", COVERAGE))
 coverage_upper_boundary <- as.numeric(sub(".*_", "", COVERAGE))
 coverage_sample_range <- seq(coverage_lower_boundary, coverage_upper_boundary, by = 1)
 set.seed(SEED)
 coverage_sample <- sample(coverage_sample_range, size = length(transcriptome), replace = T)
 per_transcript_read_N <- round(coverage_sample * width(transcriptome) / READ_LENGTH)
 expr_transcriptome_path <- TRANSCRIPTOME_PATH 
} else if (SIM_TYPE == "realistic") {
 exp_prop_lower_boundary <- as.numeric(sub("_.*", "", EXPRESSED_PROP_RANGE))
 exp_prop_upper_boundary <- as.numeric(sub(".*_", "", EXPRESSED_PROP_RANGE))
 exp_prop_sample_range <- seq(exp_prop_lower_boundary, exp_prop_upper_boundary, by = 0.01)
 set.seed(SEED)
 expressed_tr_prop_sample <- sample(exp_prop_sample_range, size = 1)
 expressed_tr_N <- ceiling(length(transcriptome) * expressed_tr_prop_sample)
 set.seed(SEED)
 expressed_transcripts <- sample(1 : length(transcriptome), size = expressed_tr_N)
 transcriptome <- transcriptome[expressed_transcripts]
 frag_N_library_total <- as.numeric(COVERAGE) * sum(width(transcriptome)) / (2 * READ_LENGTH)
 set.seed(SEED)
 per_transcript_abundance <- rlnorm(length(transcriptome), meanlog = 1, sdlog = 2)
 per_transcript_abundance <- per_transcript_abundance * width(transcriptome)
 per_transcript_read_N <- round(2 * frag_N_library_total * per_transcript_abundance / sum(per_transcript_abundance))
 expr_transcriptome_path <- paste(SIM_TYPE, "sim", OUT_PREF, READ_LENGTH, COVERAGE, MISMATCH_RATE, expressed_tr_prop_sample, SEED, "tmp_transcriptome.fasta", sep = "_")
 writeXStringSet(transcriptome, expr_transcriptome_path, format = "fasta")
 coverage_dt <- data.table(transcript = names(transcriptome), coverage = per_transcript_read_N * READ_LENGTH / width(transcriptome))
 write.table(coverage_dt, file = paste(paste(SIM_TYPE, "sim", OUT_PREF, READ_LENGTH, COVERAGE, MISMATCH_RATE, expressed_tr_prop_sample, SEED, "coverage_table.tsv", sep = "_")), sep = "\t", row.names = F, col.names = F, quote = F)
}
per_transcript_read_N[per_transcript_read_N == 0] <- 1

#Defining the fold change matrix
fold_change_matrix <- matrix(rep(1, length(transcriptome)))

#Simulating reads. Writing reads to files
set.seed(SEED)
simulate_experiment(fasta = expr_transcriptome_path, num_reps = 1, fold_changes = fold_change_matrix, paired = T, reads_per_transcript = per_transcript_read_N / 2, readlen = READ_LENGTH, error_rate = MISMATCH_RATE, distr = "empirical", bias = "rnaf", error_model = "uniform", seed = SEED, gzip = F, outdir = paste(SIM_TYPE, "sim", OUT_PREF, READ_LENGTH, COVERAGE, MISMATCH_RATE, expressed_tr_prop_sample, SEED, sep = "_"))