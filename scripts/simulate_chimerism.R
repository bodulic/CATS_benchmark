#!/usr/bin/env Rscript
#Chimerism simulation script
#Loading the required package
suppressPackageStartupMessages(library(Biostrings))
options(scipen = 999)

#Defining a sampling function that intuitively handles n = 1 scenarios
resample <- function(x, ...) x[sample.int(length(x), ...)]

#Defining external arguments
ext_args <- commandArgs(trailingOnly = T)
TRANSCRIPTOME_PATH <- ext_args[1]
MUT_FACTOR <- as.numeric(ext_args[2])
SEED <- as.numeric(ext_args[3])
MERGE_DELETION_LEN_PROP_MAX <- as.numeric(ext_args[4])
OUT_PREF <- ext_args[5]

#Importing transcriptome assembly from file
transcriptome <- readDNAStringSet(TRANSCRIPTOME_PATH, format = "fasta")
names(transcriptome) <- sub(" .*", "", names(transcriptome))

#Calculating the number of transcripts to merge
tr_to_merge_N <- round(2 * (MUT_FACTOR * length(transcriptome)) / (1 + MUT_FACTOR), 0)

#Defining transcripts to merge
set.seed(SEED)
tr_to_merge <- transcriptome[sample(1 : length(transcriptome), size = tr_to_merge_N)]

if(length(tr_to_merge) %% 2 == 1) {
 tr_to_merge <- tr_to_merge[-length(tr_to_merge)]
}

#Introducing chimerism (merging transcripts)
tr_to_merge_names <- names(tr_to_merge)
merged_tr <- DNAStringSet()
n <- 1

while (length(tr_to_merge) != 0) {
 set.seed(SEED)
 first_tr_index <- sample(1 : length(tr_to_merge), size = 1)
 first_tr <- tr_to_merge[first_tr_index]
 first_tr_name <- names(first_tr)
 first_tr_len <- width(first_tr)
 first_tr_end_boundary <- floor((1 - MERGE_DELETION_LEN_PROP_MAX) * first_tr_len) + 1
 set.seed(SEED + n)
 first_tr_end_sample_range <- seq(first_tr_end_boundary, first_tr_len, by = 1)
 first_tr_end_sample <- sample(first_tr_end_sample_range, size = 1) 
 first_tr <- unlist(first_tr)[1 : first_tr_end_sample]
  
 set.seed(SEED)
 second_tr_index <- resample((1 : length(tr_to_merge))[-first_tr_index], size = 1)
 second_tr <- tr_to_merge[second_tr_index]
 second_tr_name <- names(second_tr)
 second_tr_len <- width(second_tr)
 second_tr_start_boundary <- ceiling(MERGE_DELETION_LEN_PROP_MAX * second_tr_len)
 set.seed(SEED + n + 1000000)
 second_tr_start_sample_range <- seq(1, second_tr_start_boundary, by = 1)
 second_tr_start_sample <- sample(second_tr_start_sample_range, size = 1)
 second_tr <- unlist(second_tr)[second_tr_start_sample : second_tr_len]
 
 tr_to_merge <- tr_to_merge[c(-first_tr_index, -second_tr_index)]
 merged_tr[[n]] <- paste0(first_tr, second_tr)
 names(merged_tr)[n] <- paste(first_tr_name, second_tr_name, sep = "%")
  
 n <- n + 1
}

#Merging native and mutated transcripts
transcriptome <- transcriptome[names(transcriptome) %in% tr_to_merge_names == F]
names(merged_tr) <- paste(names(merged_tr), "chimeric_mut", sep = "_")
transcriptome <- c(transcriptome, merged_tr)

#Writing mutated transcriptome assembly to file
writeXStringSet(transcriptome, filepath = paste(OUT_PREF, "chimerism", MUT_FACTOR, SEED, sep = "_"), format = "fasta", append = F, compress = F, compression_level = NA)