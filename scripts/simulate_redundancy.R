#!/usr/bin/env Rscript
#Redundancy simulation script
#Loading the required package
suppressPackageStartupMessages(library(Biostrings))
options(scipen = 999)

#Defining external arguments
ext_args <- commandArgs(trailingOnly = T)
TRANSCRIPTOME_PATH <- ext_args[1]
SEED <- as.numeric(ext_args[2])
MUT_FACTOR <- as.numeric(ext_args[3])
RED_FRAC_LOWER_BOUNDARY <- as.numeric(ext_args[4])
RED_FRAC_UPPER_BOUNDARY <- as.numeric(ext_args[5])
OUT_PREF <- ext_args[6]

#Importing transcriptome assembly from file
transcriptome <- readDNAStringSet(TRANSCRIPTOME_PATH, format = "fasta")
names(transcriptome) <- sub(" .*", "", names(transcriptome))

#Defining transcripts to duplicate
set.seed(SEED) 
tr_to_duplicate <- transcriptome[sample(1 : length(transcriptome), size = MUT_FACTOR * length(transcriptome))]

#Introducing redundancy (duplication)
duplicated_tr_seg <- DNAStringSet()  
for (i in 1 : length(tr_to_duplicate)) {
 set.seed(SEED + i) 
 red_frac_sample_range <- seq(from = RED_FRAC_LOWER_BOUNDARY, to = RED_FRAC_UPPER_BOUNDARY, by = 0.01)
 red_frac_sample <- sample(red_frac_sample_range, size = 1)
 transcript <- tr_to_duplicate[[i]] 
 red_frac <- round(red_frac_sample * length(transcript), 0)
 transcript <- transcript[1 : red_frac]
 duplicated_tr_seg[[i]] <- transcript
}

#Merging native and mutated transcripts
transcriptome <- transcriptome[names(transcriptome) %in% names(tr_to_duplicate) == F]
names_tr_to_duplicate_tmp <- names(tr_to_duplicate)
names(tr_to_duplicate) <- paste(names(tr_to_duplicate), "redundancy_mut_native", sep = "_")
names(duplicated_tr_seg) <- paste(names_tr_to_duplicate_tmp, "redundancy_mut_added", sep = "_")
transcriptome <- c(transcriptome, tr_to_duplicate, duplicated_tr_seg)

#Writing mutated transcriptome assembly to file
writeXStringSet(transcriptome, filepath = paste(OUT_PREF, "redundancy", MUT_FACTOR, SEED, sep = "_"), format = "fasta", append = F, compress = F, compression_level = NA)