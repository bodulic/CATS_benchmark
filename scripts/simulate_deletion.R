#!/usr/bin/env Rscript
#Deletion simulation script
#Loading the required package
suppressPackageStartupMessages(library(Biostrings))
options(scipen = 999)

#Defining external arguments
ext_args <- commandArgs(trailingOnly = T)
TRANSCRIPTOME_PATH <- ext_args[1]
SEED <- as.numeric(ext_args[2])
MUT_FACTOR <- as.numeric(ext_args[3])
REGION_N <- as.numeric(ext_args[4])
OUT_PREF <- ext_args[5]

#Importing transcriptome assembly from file
transcriptome <- readDNAStringSet(TRANSCRIPTOME_PATH, format = "fasta")
names(transcriptome) <- sub(" .*", "", names(transcriptome))

#Defining transcripts to mutate
set.seed(SEED) 
tr_to_mutate <- transcriptome[sample(1 : length(transcriptome), size = MUT_FACTOR * length(transcriptome))]

#Introducing deletions
for (i in 1 : length(tr_to_mutate)) {
 transcript <- tr_to_mutate[[i]]

#Determining the number and length of introduced deletions
 set.seed(SEED + i)  
 region_N <- sample(1 : REGION_N, size = 1)
 region_length <- round((MUT_FACTOR * length(transcript)) / region_N, 0)
 
#Determining deletion coordinates while excluding already introduced deletions (in case of multiple introductions)
 exclude_range <- c()
 for (j in 1 : region_N) {
  if (j > 1) {
   exclude_range <- c(exclude_range, seq(from = start_sample - region_length + 1, to = start_sample + region_length - 1, by = 1))
  }
  start_sample_range <- 1 : (length(transcript) - region_length + 1)
  set.seed(SEED + i + j)  
  start_sample <- sample(start_sample_range[start_sample_range %in% exclude_range == F], size = 1)
  end_sample <- start_sample + region_length - 1
  transcript <- transcript[-(start_sample : end_sample)]
  }
 tr_to_mutate[[i]] <- transcript
}

#Merging native and mutated transcripts
transcriptome <- transcriptome[names(transcriptome) %in% names(tr_to_mutate) == F]
names(tr_to_mutate) <- paste(names(tr_to_mutate), "deletion_mut", sep = "_")
transcriptome <- c(transcriptome, tr_to_mutate)

#Writing mutated transcriptome assembly to file
writeXStringSet(transcriptome, filepath = paste(OUT_PREF, "deletion", MUT_FACTOR, SEED, sep = "_"), format = "fasta", append = F, compress = F, compression_level = NA)