# CATS Benchmark

This repository contains the scripts used in the analysis presented in the paper: “Comprehensive Transcriptome Quality Assessment Using CATS: Reference-free and Reference-based Approaches” (Bodulić and Vlahoviček, 2025). More information can be found in the [preprint](https://www.biorxiv.org/content/10.1101/2025.07.22.666112v1).

# General Information

To run the full pipeline, clone this repository and execute the master script (CATS_benchmark.bash) within the repository environment:

```bash
./CATS_benchmark.bash
```
Make sure all required dependencies are added to the `PATH` environment variable. All scripts in the scripts directory should also be added to `PATH` before running the script. An active internet connection is required.

# Dependencies

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
| TransRate        | 1.04        |
| seqkit           | 2.3.0       |
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
| cowplot (R)      | 1.1.3       |
| ggbiplot (R)     | 0.6.2       |

Tools denoted with (R) correspond to R packages. Version information for dependencies specific to CATS-rf and CATS-rb can be found in the README files of their respective repositories.

The pipeline was run on a high-performance computing cluster using SUSE Linux Enterprise High Performance Computing 15 SP4 (SLE_HPC 15-SP4) with Linux kernel 5.14.21. The system used a GNU Bash shell version 4.4.23, and ran on a node equipped with Intel(R) Xeon(R) Gold 6248 CPUs @ 2.50GHz. 

# Script overview

This repository contains the following scripts:

| **Script**                          | **Purpose**                                                 |        
|-------------------------------------|-------------------------------------------------------------|
| `CATS_benchmark.bash`               | Master script for the analysis                              |
| `simulate_reads.R`                  | Performs RNA-seq library simulation                         |
| `get_f_scores_from_crb_table.R`     | Calculated transcript F-scores from CRB-BLAST results       |
| `simulate_insertion.R`              | Simulates insertions in transcript sequences                |
| `simulate_mismatch.R`               | Simulates mismatches in transcript sequences                |
| `simulate_deletion.R`               | Simulates deletions in transcript sequences                 |
| `simulate_redundancy.R`             | Simulates redundancy in transcript sequences                |
| `simulate_fragmentation.R`          | Simulates fragmentation in transcript sequences             |
| `simulate_chimerism.R`              | Simulates chimerism in transcript sequences                 |
| `generate_grid_search_parameters.R` | Performs random grid search for CATS-rf parameters          |
| `filter_by_blat.R`                  | Filters public transcripts by F-scores obtained by blat     |
| `merge_cats_rf_results.R`           | Merges all CATS-rf results for Figure 2                     |
| `generate_fig2_elements.R`          | Generates all Figure 2 panels and Supplementary figures 1-3 |
| `generate_fig3_elements.R`          | Generates Figure 3 panels                                   |
| `merge_cats_rb_results.R`           | Merges all CATS-rb results for Figure 5                     |
| `generate_fig5_elements.R`          | Generates Figure 5 panels and Supplementary figures 4-9     |
| `generate_supp_fig10.R`             | Generates Supplementary Figure 10                           |
| `generate_supp_fig11.R`             | Generates Supplementary Figure 11                           |
| `generate_supp_fig12.R`             | Generates Supplementary Figure 12                           |

The repository also contains configfiles for SOAPdenovo-Trans runs for simulated (`soap_configfile_simulated`) and public (`soap_configfile_public`) libraries.

For convenience, tables containing benchmark results used in figure generation are also provided:

| **Table**                                                                            | **Figure**                        | **Download link**                                                                                                                                       | 
|--------------------------------------------------------------------------------------|-----------------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------|
| `merged_controlled_simulated_transcript_scores_for_figure2.tsv`                      | Figure 2 (a–f), Supp. Figures 1-2 | [Download](http://hex.bioinfo.hr/~kbodulic/CATS_benchmark_tables/merged_controlled_simulated_transcript_scores_for_figure2.tsv.gz)                      |
| `merged_realistically_simulated_transcript_scores_for_figure2.tsv`                   | Figure 2 (g-i)                    | [Download](http://hex.bioinfo.hr/~kbodulic/CATS_benchmark_tables/merged_realistically_simulated_transcript_scores_for_figure2.tsv.gz)                   |
| `merged_public_transcript_scores_for_figure2.tsv`                                    | Figure 2j, Supp. Figure 3         | [Download](http://hex.bioinfo.hr/~kbodulic/CATS_benchmark_tables/merged_public_transcript_scores_for_figure2.tsv.gz)                                    |
| `mutation_analysis_CATS_rf_transcript_scores_for_figure3.tsv`                        | Figure 3                          | [Download](http://hex.bioinfo.hr/~kbodulic/CATS_benchmark_tables/mutation_analysis_CATS_rf_transcript_scores_for_figure3.tsv.gz)                        |
| `merged_cats_rb_controlled_simulated_assembly_scores_for_figure5.tsv`                | Figure 5 (a–c,g,h)                | [Download](http://hex.bioinfo.hr/~kbodulic/CATS_benchmark_tables/merged_cats_rb_controlled_simulated_assembly_scores_for_figure5.tsv.gz)                |
| `merged_cats_rb_realistically_simulated_assembly_scores_for_figure5.tsv`             | Figure 5d                         | [Download](http://hex.bioinfo.hr/~kbodulic/CATS_benchmark_tables/merged_cats_rb_realistically_simulated_assembly_scores_for_figure5.tsv.gz)             |
| `merged_cats_rb_realistically_simulated_assembly_scores_per_library_for_figure5.tsv` | Figure 5e                         | [Download](http://hex.bioinfo.hr/~kbodulic/CATS_benchmark_tables/merged_cats_rb_realistically_simulated_assembly_scores_per_library_for_figure5.tsv.gz) |
| `merged_cats_rb_public_assembly_scores_for_figure5.tsv`                              | Figure 5f                         | [Download](http://hex.bioinfo.hr/~kbodulic/CATS_benchmark_tables/merged_cats_rb_public_assembly_scores_for_figure5.tsv.gz)                              |
| `controlled_simulated_assembly_cats_rf_f_scores_for_figure5.tsv`                     | Figure 5 (g,h)                    | [Download](http://hex.bioinfo.hr/~kbodulic/CATS_benchmark_tables/controlled_simulated_assembly_cats_rf_f_scores_for_figure5.tsv.gz)                     |
| `merged_cats_rb_simulated_mutated_assembly_scores_for_figure5.tsv`                   | Figure 5i                         | [Download](http://hex.bioinfo.hr/~kbodulic/CATS_benchmark_tables/merged_cats_rb_simulated_mutated_assembly_scores_for_figure5.tsv.gz)                   |
| `chimeric_ref_transcriptome_names_for_supp_figure10.tsv`                             | Supp. Figure 10                   | [Download](http://hex.bioinfo.hr/~kbodulic/CATS_benchmark_tables/chimeric_ref_transcriptome_names_for_supp_figure10.tsv.gz)                             |
| `str_inconsistent_transcripts_for_supp_figure10.tsv`                                 | Supp. Figure 10                   | [Download](http://hex.bioinfo.hr/~kbodulic/CATS_benchmark_tables/str_inconsistent_transcripts_for_supp_figure10.tsv.gz)                                 |
| `cats_rf_transcript_scores_for_supp_fig11.tsv`                                       | Supp. Figure 11                   | [Download](http://hex.bioinfo.hr/~kbodulic/CATS_benchmark_tables/cats_rf_transcript_scores_for_supp_fig11.tsv.gz)                                       |
| `transcript_f_scores_supp_fig11.tsv                `                                 | Supp. Figure 11                   | [Download](http://hex.bioinfo.hr/~kbodulic/CATS_benchmark_tables/transcript_f_scores_supp_fig11.tsv.gz)                                                 |
| `coverage_table_for_supp_figure12.tsv`                                               | Supp. Figure 12                   | [Download](http://hex.bioinfo.hr/~kbodulic/CATS_benchmark_tables/coverage_table_for_supp_figure12.tsv.gz)                                               |
| `ref_tr_size_for_supp_figure12.tsv`                                                  | Supp. Figure 12                   | [Download](http://hex.bioinfo.hr/~kbodulic/CATS_benchmark_tables/ref_tr_size_for_supp_figure12.tsv.gz)                                                  |

Tables should be unzipped before being directly supplied to the corresponding R scripts.

If the files are not directly accessible via the provided links, they are available in the following [folder](http://hex.bioinfo.hr/~kbodulic/CATS_benchmark_tables).

The tables have also been deposited to Zenodo (10.5281/zenodo.18020437)
