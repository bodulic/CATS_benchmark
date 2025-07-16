#!/usr/bin/env Rscript
#Read simulation script
#Loading the required packages
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(polyester))
options(scipen = 999)

#Defining external arguments
ext_args <- commandArgs(trailingOnly = T)
TRANSCRIPTOME_PATH <- ext_args[1]
COVERAGE_RANGE <- ext_args[2]
SEED <- as.numeric(ext_args[3])
READ_LENGTH <- as.numeric(ext_args[4])
MISMATCH_RATE <- as.numeric(ext_args[5])
OUT_PREF <- ext_args[6]

#Importing transcriptome assembly from file
transcriptome <- readDNAStringSet(TRANSCRIPTOME_PATH, format = "fasta")

#Defining coverage range
coverage_lower_boundary <- as.numeric(sub("_.*", "", COVERAGE_RANGE))
coverage_upper_boundary <- as.numeric(sub(".*_", "", COVERAGE_RANGE))
coverage_sample_range <- seq(coverage_lower_boundary, coverage_upper_boundary, by = 1)
  
#Calculating the number of simulated reads per transcript (function of coverage, transcript length, and read length)
set.seed(SEED)
coverage_sample <- sample(coverage_sample_range, size = 1)
per_transcript_read_N <- round(coverage_sample * width(transcriptome) / READ_LENGTH)
per_transcript_read_N[per_transcript_read_N == 0] <- 1

#Defining the fold change matrix
fold_change_matrix <- matrix(rep(1, length(transcriptome)))

#Simulating reads. Writing reads to files
set.seed(SEED)
simulate_experiment(fasta = TRANSCRIPTOME_PATH, num_reps = 1, fold_changes = fold_change_matrix, paired = T, reads_per_transcript = per_transcript_read_N, readlen = READ_LENGTH, error_rate = MISMATCH_RATE, distr = "empirical", bias = "rnaf", error_model = "uniform", seed = SEED, gzip = F, outdir = paste("sim", OUT_PREF, READ_LENGTH, coverage_lower_boundary, coverage_upper_boundary, MISMATCH_RATE, sep = "_"))