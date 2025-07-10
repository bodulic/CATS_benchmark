#!/usr/bin/env Rscript
#Mismatch simulation script
#Loading the required package
suppressPackageStartupMessages(library(Biostrings))
options(scipen = 999)

#Defining external arguments
ext_args <- commandArgs(trailingOnly = T)
TRANSCRIPTOME_PATH <- ext_args[1]
SAME_BASE_PROB <- as.numeric(ext_args[2])
SEED <- as.numeric(ext_args[3])
MUT_FACTOR <- as.numeric(ext_args[4])
REGION_N <- as.numeric(ext_args[5])
OUT_PREF <- ext_args[6]

#Importing transcriptome assembly from file
transcriptome <- readDNAStringSet(TRANSCRIPTOME_PATH, format = "fasta")
names(transcriptome) <- sub(" .*", "", names(transcriptome))

#Defining substitution probabilities
prob_matrix <- matrix(c(SAME_BASE_PROB, (1 - SAME_BASE_PROB) / 3, (1 - SAME_BASE_PROB) / 3, (1 - SAME_BASE_PROB) / 3,
                        (1 - SAME_BASE_PROB) / 3, SAME_BASE_PROB, (1 - SAME_BASE_PROB) / 3, (1 - SAME_BASE_PROB) / 3,
                        (1 - SAME_BASE_PROB) / 3, (1 - SAME_BASE_PROB) / 3, SAME_BASE_PROB, (1 - SAME_BASE_PROB) / 3,
                        (1 - SAME_BASE_PROB) / 3, (1 - SAME_BASE_PROB) / 3, (1 - SAME_BASE_PROB) / 3, SAME_BASE_PROB), nrow = 4, byrow = T)
nucleotides <- c("A", "T", "G", "C")

#Removing transcripts with unknown bases (Ns)
orig_transcriptome_size <- length(transcriptome)
transcriptome <- transcriptome[grepl("N", transcriptome) == F]

#Defining transcripts to mutate
set.seed(SEED)
tr_to_mutate <- transcriptome[sample(1 : length(transcriptome), size = MUT_FACTOR * orig_transcriptome_size)]

#Introducing mismatches
for (i in 1 : length(tr_to_mutate)) {
 transcript <- tr_to_mutate[[i]] 

#Determining the number and length of introduced mismatching regions
 set.seed(SEED + i)  
 region_N <- sample(1 : REGION_N, size = 1)
 region_length <- round((MUT_FACTOR * length(transcript)) / region_N, 0)

#Determining mismatching region coordinates while excluding already introduced mismatching regions (in case of multiple introductions)
 exclude_range <- c()
 for (j in 1 : region_N) {
  if (j > 1) {
   exclude_range <- c(exclude_range, seq(from = start_sample - region_length + 1, to = start_sample + region_length - 1, by = 1))
  }
  start_sample_range <- 1 : (length(transcript) - region_length + 1)
  set.seed(SEED + i + j)  
  start_sample <- sample(start_sample_range[start_sample_range %in% exclude_range == F], size = 1)
  end_sample <- start_sample + region_length - 1
  region_to_mutate <- unlist(strsplit(as.character(transcript[start_sample : end_sample]), split = ""))
  for (k in 1 : length(region_to_mutate)) {
   base_index <- match(region_to_mutate[k], nucleotides) 
   set.seed(SEED + i + j + k)  
   transcript[start_sample : end_sample][k] <- sample(nucleotides, size = 1, prob = prob_matrix[base_index, ])
  }
 }
 tr_to_mutate[[i]] <- transcript
}

#Merging native and mutated transcripts
transcriptome <- transcriptome[names(transcriptome) %in% names(tr_to_mutate) == F]
names(tr_to_mutate) <- paste(names(tr_to_mutate), "mismatch_mut", sep = "_")
transcriptome <- c(transcriptome, tr_to_mutate)

#Writing mutated transcriptome assembly to file
writeXStringSet(transcriptome, filepath = paste(OUT_PREF, "mismatch", MUT_FACTOR, SEED, sep = "_"), format = "fasta", append = F, compress = F, compression_level = NA)