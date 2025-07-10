#!/usr/bin/env Rscript
#Fragmentation simulation script
#Loading the required package
suppressPackageStartupMessages(library(Biostrings))
options(scipen = 999)

#Defining external arguments
ext_args <- commandArgs(trailingOnly = T)
TRANSCRIPTOME_PATH <- ext_args[1]
SEED <- as.numeric(ext_args[2])
MUT_FACTOR <- as.numeric(ext_args[3])
FRAG_POINT_LOWER_BOUNDARY <- as.numeric(ext_args[4])
FRAG_POINT_UPPER_BOUNDARY <- as.numeric(ext_args[5])
OUT_PREF <- ext_args[6]

#Importing transcriptome assembly from file
transcriptome <- readDNAStringSet(TRANSCRIPTOME_PATH, format = "fasta")
names(transcriptome) <- sub(" .*", "", names(transcriptome))

#Defining transcripts to fragment
set.seed(SEED) 
tr_to_fragment <- transcriptome[sample(1 : length(transcriptome), size = MUT_FACTOR * length(transcriptome))]

#Introducing fragmentation
first_tr_frag <- DNAStringSet()  
second_tr_frag <- DNAStringSet()  
for (i in 1 : length(tr_to_fragment)) {
 set.seed(SEED + i)  
 frag_point_sample_range <- seq(from = FRAG_POINT_LOWER_BOUNDARY, to = FRAG_POINT_UPPER_BOUNDARY, by = 0.01)
 frag_point_sample <- sample(frag_point_sample_range, size = 1)
 transcript <- tr_to_fragment[[i]]
 frag_point <- round(frag_point_sample * length(transcript), 0)
 first_tr_frag[[i]] <- transcript[1 : frag_point]
 second_tr_frag[[i]] <- transcript[(frag_point + 1) : length(transcript)]
}

#Merging native and mutated transcripts
transcriptome <- transcriptome[names(transcriptome) %in% names(tr_to_fragment) == F]
names(first_tr_frag) <- paste(names(tr_to_fragment), "fragmentation_mut1", sep = "_")
names(second_tr_frag) <- paste(names(tr_to_fragment), "fragmentation_mut2", sep = "_")
transcriptome <- c(transcriptome, first_tr_frag, second_tr_frag)

#Writing mutated transcriptome assembly to file
writeXStringSet(transcriptome, filepath = paste(OUT_PREF, "fragmentation", MUT_FACTOR, SEED, sep = "_"), format = "fasta", append = F, compress = F, compression_level = NA)