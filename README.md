# CATS Bemchmrk

This repository contains the scripts used in the analyses presented in the paper: “Comprehensive Transcriptome Quality Assessment Using CATS: Reference-free and Reference-based Approaches” (Bodulić & Vlahoviček, 2025). More information can be found in the [preprint](x).

# General Information

To run the full pipeline, clone this repository and execute the master script (CATS_benchmark.bash) within the repository environment:

```bash
./CATS_benchmark.bash
```

Make sure all required dependencies are available in your environment before running the script. Active internet connection is required.

## Dependencies

The following dependencies were used in the pipeline:

| **Dependency**   | **Version** |                                                              
|------------------|-------------|
| R                | 4.4.2       |
| Biostrings (R)   | 2.27.1      |
| polyester (R)    | 1.32.0      |
| SPAdes           | 3.15.4      |
| Trinity          | 2.15.1      |
| seqtk            | 1.3-r106    |
| IDBA-tran        | 1.1.1       |
| SOAPdenovo-Trans | 1.04        |
| CATS-rf          | 1.0.0       |
| RSEM-EVAL        | 1.11        |
| Ruby             | 2.6.6p146   |
| TramsRate        | 1.04        |
| BLASTn           | 2.2.29      |
| CRB-BLAST        | 1.0.0       |
| data.table (R)   | 1.16.4      |
| CATS-rb          | 1.0.0       |
| SRA toolkit      | 3.2.0       |
| Python           | 3.10.6      |
| cutadapt         | 4.6         |
| TrimGalore       | 0.6.1       |
| kallisto         | 0.50.1      |
| pblat            | 2.5.1       |
| ggplot2 (R)      | 3.5.1       | 
| ggsci (R)        | 3.2.0       |
| ggcorrplot (R)   | 0.1.4.1     |
| ggbiplot (R)     | 0.6.2       |
| cowplot (R)      | 1.1.3       |

The stated dependencies should be included in `PATH` variable before running the script. Tools denoted with (R) correspond to R packages. Version information for dependencies specific to CATS-rf and CATS-rb can be found in the README files of their respective repositories.

The pipeline was run on a high-performance computing cluster using SUSE Linux Enterprise High Performance Computing 15 SP4 (SLE_HPC 15-SP4) with Linux kernel 5.14.21. The system used a GNU Bash shell version 4.4.23, and ran on a node equipped with Intel(R) Xeon(R) Gold 6248 CPUs @ 2.50GHz. 

# Script overview

| **Script**                      | **Purpose**                                            |        
|---------------------------------|--------------------------------------------------------|
| `CATS_benchmark.bash`           | Master script for the analysis                         |
| `simulate_reads.R`              | Performs RNA-seq library simulation                    |
| `get_f_scores_from_crb_table.R` | Calculated transcript F-scores from CRB-BLAST results  |
| `simulate_insertion.R`          | Simulates insertions in transcript sequences           |
| `simulate_mismatch.R`           | Simulates mismatches in transcript sequences           |
| `simulate_deletion.R`           | Simulates deletions in transcript sequences            |
| `simulate_redundancy.R`         | Simulates redundancy in transcript sequences           |
| `simulate_fragmentation.R`      | Simulates fragmentation in transcript sequences        |
| `simulate_chimerism.R`          | Simulates chimerism in transcript sequences            |
| `filter_by_blat`                | Filters public transcripts by blat F-scores            |
| `merge_cats_rf_results.R`       | Merges all CATS-rf results for Figure 2                |
| `generate_fig2_elements.R`      | Generates all Figure 2 subplots and Ext. data figure 1 |
| `generate_fig2_elements.R`      | Generates Figure 3 subplots                            |
| `merge_cats_rb_results.R`       | Merges all CATS-rb results for Figure 5                |
| `generate_fig5_elements.R`      | Generates Figure 5 subplots and Ext. data figures 2-6  |
| `generate_fig_ext7.R`             | Generates Ext. data Figure 7                           |

For convinience, tables containing benchmark results used in figure generation are provided:

| **Table**                                              | **Figure**                    | **Download link**                                                                                                                     | 
|--------------------------------------------------------|-------------------------------|---------------------------------------------------------------------------------------------------------------------------------------|
| `merged_simulated_transcript_scores_for_figure2.tsv`   | Figure 2 (a–e)                | [Download](http://hex.bioinfo.hr/~kbodulic/CATS_benchmark_tables/merged_simulated_transcript_scores_for_figure2.tsv.gz)               |
| `merged_simulated_assembly_scores_for_figure2.tsv`     | Figure 2f                     | [Download](http://hex.bioinfo.hr/~kbodulic/CATS_benchmark_tables/merged_simulated_assembly_scores_for_figure2.tsv.gz)                 |
| `merged_public_transcript_scores_for_figure2.tsv`      | Figure 2g, Ext data. figure 1 | [Download](http://hex.bioinfo.hr/~kbodulic/CATS_benchmark_tables/merged_public_transcript_scores_for_figure2.tsv.gz)                  |
| `mutation_analysis_CATS_rf_transcript_scores.tsv`      | Figure 3                      | [Download](http://hex.bioinfo.hr/~kbodulic/CATS_benchmark_tables/mutation_analysis_CATS_rf_transcript_scores_for_figure3.tsv.gz)      |
| `merged_cats_rb_simulated_assembly_scores.tsv`         | Figure 5 (a–c,e,f)            | [Download](http://hex.bioinfo.hr/~kbodulic/CATS_benchmark_tables/merged_cats_rb_simulated_assembly_scores_for_figure5.tsv.gz)         |
| `merged_cats_rb_public_assembly_scores.tsv`            | Figure 5d                     | [Download](http://hex.bioinfo.hr/~kbodulic/CATS_benchmark_tables/merged_cats_rb_public_assembly_scores_for_figure5.tsv.gz)            |
| `simulated_assembly_cats_rf_f_scores.tsv`              | Figure 5 (e,f)                | [Download](http://hex.bioinfo.hr/~kbodulic/CATS_benchmark_tables/simulated_assembly_cats_rf_f_scores_for_figure5.tsv.gz)              |
| `merged_cats_rb_simulated_mutated_scores.tsv`          | Figure 5g                     | [Download](http://hex.bioinfo.hr/~kbodulic/CATS_benchmark_tables/merged_cats_rb_simulated_mutated_assembly_scores_for_figure5.tsv.gz) |
| `chimeric_ref_transcriptome_names.tsv`                 | Ext. data figure 7            | [Download](http://hex.bioinfo.hr/~kbodulic/CATS_benchmark_tables/chimeric_ref_transcriptome_names_for_figure_ext7.gz)                   |
| `str_inconsistent_transcripts.tsv`                     | Ext. data figure 7            | [Download](http://hex.bioinfo.hr/~kbodulic/CATS_benchmark_tables/str_inconsistent_transcripts_for_figure_ext7.gz)                       |

Tables should be unzipped before being directly supplied to the correpsonding R scripts.
