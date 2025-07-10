#!/bin/bash
#Master script for the analysis performed in the paper "Comprehensive Transcriptome Quality Assessment Using CATS: Reference-free and Reference-based Approaches" (Bodulić and Vlahoviček, 2025)

#Downloading reference protein-coding and non-coding transcripts from Ensembl
declare -A transcriptome_urls
declare -A ref_tr_outputs

#S. cerevisiae
transcriptome_urls[s_cerevisiae_cdna]="https://ftp.ensembl.org/pub/release-111/fasta/saccharomyces_cerevisiae/cdna/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz"
transcriptome_urls[s_cerevisiae_ncrna]="https://ftp.ensembl.org/pub/release-111/fasta/saccharomyces_cerevisiae/ncrna/Saccharomyces_cerevisiae.R64-1-1.ncrna.fa.gz"
ref_tr_outputs[s_cerevisiae]="s_cerevisiae_r64_1_1_transcripts.fa"

#C. elegans
transcriptome_urls[c_elegans_cdna]="https://ftp.ensembl.org/pub/release-111/fasta/caenorhabditis_elegans/cdna/Caenorhabditis_elegans.WBcel235.cdna.all.fa.gz"
transcriptome_urls[c_elegans_ncrna]="https://ftp.ensembl.org/pub/release-111/fasta/caenorhabditis_elegans/ncrna/Caenorhabditis_elegans.WBcel235.ncrna.fa.gz"
ref_tr_outputs[c_elegans]="c_elegans_wbcel235_transcripts.fa"

#D. melanogaster
transcriptome_urls[d_melanogaster_cdna]="https://ftp.ensembl.org/pub/release-111/fasta/drosophila_melanogaster/cdna/Drosophila_melanogaster.BDGP6.46.cdna.all.fa.gz"
transcriptome_urls[d_melanogaster_ncrna]="https://ftp.ensembl.org/pub/release-111/fasta/drosophila_melanogaster/ncrna/Drosophila_melanogaster.BDGP6.46.ncrna.fa.gz"
ref_tr_outputs[d_melanogaster]="d_melanogaster_bdgp6_46_transcripts.fa"

#A. thaliana
transcriptome_urls[a_thaliana_cdna]="http://ftp.ensemblgenomes.org/pub/plants/release-56/fasta/arabidopsis_thaliana/cdna/Arabidopsis_thaliana.TAIR10.cdna.all.fa.gz"
transcriptome_urls[a_thaliana_ncrna]="http://ftp.ensemblgenomes.org/pub/plants/release-56/fasta/arabidopsis_thaliana/ncrna/Arabidopsis_thaliana.TAIR10.ncrna.fa.gz"
ref_tr_outputs[a_thaliana]="a_thaliana_tair10_transcripts.fa"

#M. musculus
transcriptome_urls[m_musculus_cdna]="https://ftp.ensembl.org/pub/release-111/fasta/mus_musculus/cdna/Mus_musculus.GRCm39.cdna.all.fa.gz"
transcriptome_urls[m_musculus_ncrna]="https://ftp.ensembl.org/pub/release-111/fasta/mus_musculus/ncrna/Mus_musculus.GRCm39.ncrna.fa.gz"
ref_tr_outputs[m_musculus]="m_musculus_grcm39_transcripts.fa"

#H. sapiens
transcriptome_urls[h_sapiens_cdna]="https://ftp.ensembl.org/pub/release-111/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz"
transcriptome_urls[h_sapiens_ncrna]="https://ftp.ensembl.org/pub/release-111/fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz"
ref_tr_outputs[h_sapiens]="h_sapiens_grhc38_transcripts.fa"

mkdir ref_transcriptomes && cd ref_transcriptomes
for key in s_cerevisiae c_elegans d_melanogaster a_thaliana m_musculus h_sapiens
do
 cdna_file="${key}_cdna.fa.gz"
 ncrna_file="${key}_ncrna.fa.gz"
 wget -O "${cdna_file}" "${transcriptome_urls[${key}_cdna]}"
 wget -O "${ncrna_file}" "${transcriptome_urls[${key}_ncrna]}"
 gunzip -c "${cdna_file}" > "${ref_tr_outputs[${key}]}"
 gunzip -c "${ncrna_file}" >> "${ref_tr_outputs[${key}]}"
 rm "${cdna_file}" "${ncrna_file}"
done
cd ..

# Downloading reference genomes from Ensembl
genome_urls=(
 "https://ftp.ensembl.org/pub/release-111/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz"
 "https://ftp.ensembl.org/pub/release-111/fasta/caenorhabditis_elegans/dna/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa.gz"
 "https://ftp.ensembl.org/pub/release-111/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.46.dna.toplevel.fa.gz"
 "http://ftp.ensemblgenomes.org/pub/plants/release-56/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz"
 "https://ftp.ensembl.org/pub/release-111/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz"
 "https://ftp.ensembl.org/pub/release-111/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
)

mkdir ref_genomes && cd ref_genomes
for genome_url in "${genome_urls[@]}"
do
 output="$(basename "${genome_url}")"
 wget "${genome_url}"
 gunzip "${output}"
done
cd ..

# Downloading reference GTF annotation from Ensembl. Removing GTF records not found in the corresponding reference transcriptome
gtf_urls=(
 "https://ftp.ensembl.org/pub/release-111/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.111.gtf.gz"
 "https://ftp.ensembl.org/pub/release-111/gtf/caenorhabditis_elegans/Caenorhabditis_elegans.WBcel235.111.gtf.gz"
 "https://ftp.ensembl.org/pub/release-111/gtf/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.46.111.gtf.gz"
 "http://ftp.ensemblgenomes.org/pub/plants/release-56/gtf/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.56.gtf.gz"
 "https://ftp.ensembl.org/pub/release-111/gtf/mus_musculus/Mus_musculus.GRCm39.111.gtf.gz"
 "https://ftp.ensembl.org/pub/release-111/gtf/homo_sapiens/Homo_sapiens.GRCh38.111.gtf.gz"
)

mkdir ref_annotation && cd ref_annotation
for gtf_url in "${gtf_urls[@]}"
do
 output="$(basename "${gtf_url}")"
 wget "${gtf_url}"
 gunzip "${output}"
 case "${output}" in
  Saccharomyces_cerevisiae.R64-1-1.111.gtf.gz) ref_transcriptome="${ref_tr_outputs[s_cerevisiae]}" ;;
  Caenorhabditis_elegans.WBcel235.111.gtf.gz)  ref_transcriptome="${ref_tr_outputs[c_elegans]}" ;;
  Drosophila_melanogaster.BDGP6.46.111.gtf.gz) ref_transcriptome="${ref_tr_outputs[d_melanogaster]}" ;;
  Arabidopsis_thaliana.TAIR10.56.gtf.gz)       ref_transcriptome="${ref_tr_outputs[a_thaliana]}" ;;
  Mus_musculus.GRCm39.111.gtf.gz)              ref_transcriptome="${ref_tr_outputs[m_musculus]}" ;;
  Homo_sapiens.GRCh38.111.gtf.gz)              ref_transcriptome="${ref_tr_outputs[h_sapiens]}" ;;
 esac
 grep "^>" "../ref_transcriptomes/${ref_transcriptome}" | cut -d' ' -f1 | sed 's/^>//; s/\..*//' | sort -u > .tmp_ids && awk 'BEGIN { while ((getline < ".tmp_ids") > 0) ids[$1] = 1 } { match($0, /transcript_id "([^"]+)"/, m); tid = m[1]; sub(/\..*/, "", tid); if (tid in ids) print }' "${output%.gz}" > "adj_${output%.gz}" && rm "${output%.gz}" .tmp_ids
done
cd ..

#Simulating RNA-seq reads
coverage_levels=("1_4" "5_10" "11_20" "21_30" "31_40" "41_50" "51_100")
mismatch_rates=(0.005 0.01 0.02)

mkdir simulated_data && cd simulated_data
ref_transcriptomes="$(ls -tr ../ref_transcriptomes/*)"
for ref_transcriptome in ${ref_transcriptomes}
do
 ref_transcriptome_name="$(basename "${ref_transcriptome}")"
 for coverage in "${coverage_levels[@]}"
 do
  for mismatch in "${mismatch_rates[@]}"
  do
   echo "Simulating reads from ${ref_transcriptome_name}, Coverage: ${coverage}, Read length: 100 bp, Mismatch rate: ${mismatch}"
   simulate_reads.R "${ref_transcriptome}" "${coverage}" 12345 100 "${mismatch}" "${ref_transcriptome_name}"
  done
 done
done

#Assembling simulated libraries with rnaSPAdes, Trinity, IDBA-tran, and SOAPdenovo-Trans
for dir in sim_*
do
 cd "${dir}"

 outdir="${dir}_rnaspades_dir"
 echo "Assembling ${dir} with rnaSPAdes"
 rnaspades.py -1 sample_01_1.fasta -2 sample_01_2.fasta -t 20 -m 300 -o "${outdir}"
 mv "${outdir}/transcripts.fasta" "${dir}_RSP"

 outdir="${dir}_trinity_dir"
 cat sample_01_1.fasta | awk '{print (NR%2 == 1) ? ">" ++i : $0}' > sample_01_1_triv_header.fa
 cat sample_01_2.fasta | awk '{print (NR%2 == 1) ? ">" ++i : $0}' > sample_01_2_triv_header.fa
 echo "Assembling ${dir} with Trinity"
 Trinity --left sample_01_1_triv_header.fa --right sample_01_2_triv_header.fa --seqType fa --CPU 20 --max_memory 300G --output "${outdir}"
 mv "${outdir}.Trinity.fasta" "${dir}_TRI"

 outdir="${dir}_idba_dir"
 seqtk seq -F 'I' sample_01_1.fasta > sample_01_1.fq
 seqtk seq -F 'I' sample_01_2.fasta > sample_01_2.fq
 fq2fa --merge --filter sample_01_1.fq sample_01_2.fq reads_together.fq
 echo "Assembling ${dir} with IDBA-tran"
 idba_tran -r reads_together.fq --num_threads 20 -o "${outdir}"
 mv "${outdir}/contig.fa" "${dir}_IDB"

 echo "Assembling ${dir} with SOAPdenovo-Trans"
 SOAPdenovo-Trans-31mer all -s ../../soap_configfile_simulated -p 20 -o "${dir}_soap"
 mv "${dir}_soap.scafSeq" "${dir}_SOA"

 cd ..
done

#Running CATS-rf, TransRate, and RSEM-EVAL on simulated transcriptome assemblies
for dir in sim_*
do
 cd "${dir}"
 for transcriptome in *_RSP *_TRI *_IDB *_SOA
 do

  echo "Running CATS-rf on ${transcriptome}"
  CATS_rf -t 20 -M 16384M -o "${transcriptome}_CATS_rf" "${transcriptome}" sample_01_1.fasta sample_01_2.fasta

  mkdir "${transcriptome}_rsem_eval" && cd "${transcriptome}_rsem_eval"
  echo "Running RSEM-EVAL on ${transcriptome}"
  rsem-eval-estimate-transcript-length-distribution "../${transcriptome}" parameter_file
  rsem-eval-calculate-score --seed 100 -p 20 --paired-end --transcript-length-parameters parameter_file ../sample_01_1.fq ../sample_01_2.fq "../${transcriptome}" "${transcriptome}_rsem" 180
  cut -f 2 "${transcriptome}_rsem.score" | head -n 1 > "${transcriptome}_rsem_eval_assembly_score"
  cut -f 1,9 "${transcriptome}_rsem.score.isoforms.results" > "${transcriptome}_rsem_eval_transcript_scores"
  cd ..

  echo "Running TransRate on ${transcriptome}"
  transrate --assembly="${transcriptome}" --left=sample_01_1.fq --right=sample_01_2.fq --threads=20 --output="${transcriptome}_transrate"
  cut -d "," -f 37 "${transcriptome}_transrate/assemblies.csv" | tail -n +2 > "${transcriptome}_transrate/${transcriptome}_transrate_assembly_score"
  cut -d "," -f 1,16,15,17,18,9 "${transcriptome}"_transrate/sim*/contigs.csv > "${transcriptome}_transrate/${transcriptome}_transrate_transcript_scores"

 done
 cd ..
done

#Calculating F-scores of simulated transcriptome assemblies with CRB-BLAST
map_ref_transcriptome_info() {
 case "${1}" in
  *s_cerevisiae*)   ref_transcriptome="${ref_tr_outputs[s_cerevisiae]}" ;;
  *c_elegans*)      ref_transcriptome="${ref_tr_outputs[c_elegans]}" ;;
  *d_melanogaster*) ref_transcriptome="${ref_tr_outputs[d_melanogaster]}" ;;
  *a_thaliana*)     ref_transcriptome="${ref_tr_outputs[a_thaliana]}" ;;
  *m_musculus*)     ref_transcriptome="${ref_tr_outputs[m_musculus]}" ;;
  *h_sapiens*)      ref_transcriptome="${ref_tr_outputs[h_sapiens]}";;
  esac
}

for dir in sim_*
do
 cd "${dir}"
 map_ref_transcriptome_info "${dir}"
 for transcriptome in *_RSP *_TRI *_IDB *_SOA
 do
  mkdir "${transcriptome}_crb" && cd "${transcriptome}_crb"
  echo "Running CRB-BLAST on ${transcriptome}"
  crb-blast -q "../${transcriptome}" -t "../../../ref_transcriptomes/${ref_transcriptome}" -h 20 -o "${transcriptome}_crb"
  get_f_scores_from_crb_table.R 20 "${transcriptome}_crb" "../${transcriptome}" "${transcriptome}"
  cd ..
 done
 cd ..
done

#Simulating mutations in a subset of simulated transcriptome assemblies (all species, coverage 21-30x, mismatch rate 0.05, rnaSPAdes and Trinity assemblers). Running CATS-rf on the mutated assemblies
for dir in *100_21_30_0.005
do
 cd "${dir}"
 for transcriptome in *_RSP *_TRI
 do
  mkdir "${transcriptome}_mut_sim" && cd "${transcriptome}_mut_sim"

  mkdir insertion && cd insertion
  echo "Simulating insertions in ${transcriptome}. Performing CATS-rf"
  for mut_factor in 0.15 0.3 0.45 0.6
  do
   simulate_insertion.R "../../${transcriptome}" 1 "${mut_factor}" 1 "${transcriptome}"
   CATS_rf -t 20 -M 16384M -o "${transcriptome}_CATS_rf_mut_ins_${mut_factor}" "${transcriptome}_insertion_${mut_factor}_1" ../../sample_01_1.fasta ../../sample_01_2.fasta
  done
  cd ..

  mkdir mismatch && cd mismatch
  echo "Simulating mismatches in ${transcriptome}. Performing CATS-rf"
  for mut_factor in 0.15 0.3 0.45 0.6
  do
   simulate_mismatch.R "../../${transcriptome}" 0.75 1 "${mut_factor}" 1 "${transcriptome}"
   CATS_rf -t 20 -M 16384M -o "${transcriptome}_CATS_rf_mut_mis_${mut_factor}" "${transcriptome}_mismatch_${mut_factor}_1" ../../sample_01_1.fasta ../../sample_01_2.fasta
  done
  cd ..

  mkdir deletion
  cd deletion
  echo "Simulating deletions in ${transcriptome}. Performing CATS-rf"
  for mut_factor in 0.15 0.3 0.45 0.6
  do
   simulate_deletion.R "../../${transcriptome}" 1 "${mut_factor}" 1 "${transcriptome}"
   CATS_rf -t 20 -M 16384M -o "${transcriptome}_CATS_rf_mut_del_${mut_factor}" "${transcriptome}_deletion_${mut_factor}_1" ../../sample_01_1.fasta ../../sample_01_2.fasta
  done
  cd ..

  mkdir redundancy
  cd redundancy
  red_frac_lower_bounds=(0.2 0.4 0.6 0.8)
  red_frac_upper_bounds=(0.39 0.59 0.79 0.99)
  echo "Simulating redundancy in ${transcriptome}. Performing CATS-rf"
  for i in "${!red_frac_lower_bounds[@]}"
  do
   mut_factor="${red_frac_lower_bounds[$i]}"
   red_frac_lower_bound="${red_frac_lower_bounds[$i]}"
   red_frac_upper_bound="${red_frac_upper_bounds[$i]}"
   simulate_redundancy.R "../../${transcriptome}" 1 "${mut_factor}" "${red_frac_lower_bound}" "${red_frac_upper_bound}" "${transcriptome}"
   CATS_rf -t 20 -M 16384M -o "${transcriptome}_CATS_rf_mut_red_${mut_factor}" "${transcriptome}_redundancy_${mut_factor}_1" ../../sample_01_1.fasta ../../sample_01_2.fasta
  done
  cd ..

  mkdir fragmentation
  cd fragmentation
  echo "Simulating fragmentation in ${transcriptome}. Performing CATS-rf"
  for mut_factor in 0.2 0.4 0.6 0.8
  do
   simulate_fragmentation.R "../../${transcriptome}" 1 "${mut_factor}" 0.4 0.6 "${transcriptome}"
   CATS_rf -t 20 -M 16384M -o "${transcriptome}_CATS_rf_mut_frag_${mut_factor}" "${transcriptome}_fragmentation_${mut_factor}_1" ../../sample_01_1.fasta ../../sample_01_2.fasta
  done
  cd ..

  mkdir multiple
  cd multiple
  echo "Performing multiplicative simulation in ${transcriptome}. Performing CATS-rf"
  for seed in 10000000 60000000 110000000
  do
   for mut_factor in 0.1 0.2 0.3 0.4 0.5 0.6
   do
    red_frac_lower_bound="${mut_factor}"
    red_frac_upper_bound="$(echo "${mut_factor} + 0.09" | bc)"
    simulate_deletion.R "../../${transcriptome}" "${seed}" "${mut_factor}" 1 mult1
    simulate_mismatch.R "mult1_deletion_${mut_factor}_${seed}" 0.75 "${seed}" "${mut_factor}" 1 mult2
    simulate_insertion.R "mult2_mismatch_${mut_factor}_${seed}" "${seed}" "${mut_factor}" 1 mult3
    simulate_fragmentation.R "mult3_insertion_${mut_factor}_${seed}" "${seed}" "${mut_factor}" 0.4 0.6 mult4
    simulate_redundancy.R "mult4_fragmentation_${mut_factor}_${seed}" "${seed}" "${mut_factor}" "${red_frac_lower_bound}" "${red_frac_upper_bound}" "${transcriptome}_mult5"
    CATS_rf -t 20 -M 16384M -o "${transcriptome}_CATS_rf_mut_mult_${mut_factor}_${seed}" "${transcriptome}_mult5_redundancy_${mut_factor}_${seed}" ../../sample_01_1.fasta ../../sample_01_2.fasta
   done
  done
  cd ../../
 done
 cd ..
done

#Simulating chimerism in reference transcriptomes for CATS-rb chimerism classification analysis. Only including transcripts >= 200 bp
cd ../ref_transcriptomes
ref_transcriptomes="$(ls -tr *)"
mkdir chimerism && cd chimerism

for ref_transcriptome in ${ref_transcriptomes}
do
 echo "Simulating chimerism in ${ref_transcriptome}"
 awk 'BEGIN{RS=">"; ORS=""} NR>1{h=substr($0,1,index($0,"\n")-1); s=substr($0,index($0,"\n")+1); gsub("\n","",s); if(length(s)>=200) print ">"h"\n"s"\n"}' "../${ref_transcriptome}" > "${ref_transcriptome}_len_filt"
 simulate_chimerism.R "${ref_transcriptome}_len_filt" 0.1 1 0.03 "${ref_transcriptome}"
done
cd ..

#Generating reference genome indices with CATS_rb_index
cd ../ref_genomes
ref_genomes="$(ls -tr *)"
mkdir ref_genome_indices && cd ref_genome_indices
for ref_genome in ${ref_genomes}
do
 case "${ref_genome}" in
  Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa) max_gene_size=25000 ;;
  Caenorhabditis_elegans.WBcel235.dna.toplevel.fa)  max_gene_size=100000 ;;
  Drosophila_melanogaster.BDGP6.46.dna.toplevel.fa) max_gene_size=500000 ;;
  Arabidopsis_thaliana.TAIR10.dna.toplevel.fa)      max_gene_size=50000 ;;
  Mus_musculus.GRCm39.dna.primary_assembly.fa)      max_gene_size=3000000 ;;
  Homo_sapiens.GRCh38.dna.primary_assembly.fa)      max_gene_size=3000000 ;;
 esac
 echo "Creating ${ref_genome} genome index"
 CATS_rb_index -m "${max_gene_size}" -t 20 "../${ref_genome}"
done
cd ../../

#Mapping reference transcriptomes and simulated transcriptome assemblies to the respective reference genomes with CATS_rb_map
map_species_info() {
 case "${1}" in
  *s_cerevisiae*)
   ref_genome_index=CATS_rb_index_Saccharomyces_cerevisiae_R64-1-1_dna_toplevel.fa
   species_preset=sacccere ;;
  *c_elegans*)
   ref_genome_index=CATS_rb_index_Caenorhabditis_elegans_WBcel235_dna_toplevel.fa
   species_preset=caeneleg ;;
  *d_melanogaster*)
   ref_genome_index=CATS_rb_index_Drosophila_melanogaster_BDGP6_46_dna_toplevel.fa
   species_preset=drosmela ;;
  *a_thaliana*)
   ref_genome_index=CATS_rb_index_Arabidopsis_thaliana_TAIR10_dna_toplevel.fa
   species_preset=arabthal ;;
  *m_musculus*)
   ref_genome_index=CATS_rb_index_Mus_musculus_GRCm39_dna_primary_assembly.fa
   species_preset=mus_musc ;;
  *h_sapiens*)
   ref_genome_index=CATS_rb_index_Homo_sapiens_GRCh38_dna_primary_assembly.fa
   species_preset=homosapi ;;
 esac
}

#Mapping reference transcriptomes
cd ref_transcriptomes
ref_transcriptomes="$(ls -tr *transcripts.fa)"
mkdir reference_trans_mapped && cd reference_trans_mapped
for ref_transcriptome in ${ref_transcriptomes}
do
 map_species_info "${ref_transcriptome}"
 echo "Mapping ${ref_transcriptome} to the reference genome"
 CATS_rb_map -p "${species_preset}" -t 20 "../../ref_genomes/ref_genome_indices/${ref_genome_index}" "../${ref_transcriptome}"
done
cd ..

#Mapping simulated transcriptome assemblies
cd simulated_data
for dir in sim_*
do
 cd "${dir}"
 map_species_info "${dir}"
 for transcriptome in *_RSP *_TRI *_IDB *_SOA
 do
  echo "Mapping ${transcriptome} to the reference genome"
  CATS_rb_map -p "${species_preset}" -t 20 "../../ref_genomes/ref_genome_indices/${ref_genome_index}" "${transcriptome}"
 done

#Mapping simulated transcriptome assemblies undergoing multiplicative mutation simulation
 if [[ "${dir}" == *100_21_30_0.005 ]]
 then
  for transcriptome_mut_dir in *_RSP_mut_sim *_TRI_mut_sim
  do
   cd "${transcriptome_mut_dir}/multiple"
   for transcriptome_mut in *_mult5_*
   do
    if [[ "${transcriptome_mut}" != *CATS_rf* ]]
    then
     echo "Mapping ${transcriptome_mut} to the reference genome"
     CATS_rb_map -p "${species_preset}" -t 20 "../../../../ref_genomes/ref_genome_indices/${ref_genome_index}" "${transcriptome_mut}"
    fi
   done
   cd ../../
  done
 fi
 cd ..
done

#Comparing simulated transcriptome assemblies and assemblies undergoing multiplicative mutation simulation with CATS_rb_compare
compare_species_info() {
 case "${1}" in
  *s_cerevisiae*)
   ref_genome=Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa
   ref_transcriptome=s_cerevisiae_r64_1_1_transcripts.fa
   ref_annotation=adj_Saccharomyces_cerevisiae.R64-1-1.111.gtf
   min_exon_id_prop=0.98
   max_intron_len=1500
   min_exon_set_len=0
   min_tr_set_len=100
   max_tr_set_len=25000
   num_longest_scaff=16
   num_genomic_bins=10000 ;;
  *c_elegans*)
   ref_genome=Caenorhabditis_elegans.WBcel235.dna.toplevel.fa
   ref_transcriptome=c_elegans_wbcel235_transcripts.fa
   ref_annotation=adj_Caenorhabditis_elegans.WBcel235.111.gtf
   min_exon_id_prop=0.98
   max_intron_len=20000
   min_exon_set_len=0
   min_tr_set_len=150
   max_tr_set_len=100000
   num_longest_scaff=6
   num_genomic_bins=15000 ;;
  *d_melanogaster*)
   ref_genome=Drosophila_melanogaster.BDGP6.46.dna.toplevel.fa
   ref_transcriptome=d_melanogaster_bdgp6_46_transcripts.fa
   ref_annotation=adj_Drosophila_melanogaster.BDGP6.46.111.gtf
   min_exon_id_prop=0.98
   max_intron_len=90000
   min_exon_set_len=0
   min_tr_set_len=200
   max_tr_set_len=500000
   num_longest_scaff=4
   num_genomic_bins=15000 ;;
  *a_thaliana*)
   ref_genome=Arabidopsis_thaliana.TAIR10.dna.toplevel.fa
   ref_transcriptome=a_thaliana_tair10_transcripts.fa
   ref_annotation=adj_Arabidopsis_thaliana.TAIR10.56.gtf
   min_exon_id_prop=0.98
   max_intron_len=5000
   min_exon_set_len=0
   min_tr_set_len=100
   max_tr_set_len=50000
   num_longest_scaff=5
   num_genomic_bins=15000 ;;
  *m_musculus*)
   ref_genome=Mus_musculus.GRCm39.dna.primary_assembly.fa
   ref_transcriptome=m_musculus_grcm39_transcripts.fa
   ref_annotation=adj_Mus_musculus.GRCm39.111.gtf
   min_exon_id_prop=0.998
   max_intron_len=600000
   min_exon_set_len=125
   min_tr_set_len=450
   max_tr_set_len=3000000
   num_longest_scaff=21
   num_genomic_bins=40000 ;;
  *h_sapiens*)
   ref_genome=Homo_sapiens.GRCh38.dna.primary_assembly.fa
   ref_transcriptome=h_sapiens_grhc38_transcripts.fa
   ref_annotation=adj_Homo_sapiens.GRCh38.111.gtf
   min_exon_id_prop=0.998
   max_intron_len=400000
   min_exon_set_len=100
   min_tr_set_len=300
   max_tr_set_len=3000000
   num_longest_scaff=24
   num_genomic_bins=40000 ;;
 esac
}

#Running CATS_rb_compare on the complete dataset, five assembly subsets with reducing maximum coverage and simulated assemblies undergoing multiplicative mutation simulation
mkdir CATS_rb_compare_results && cd CATS_rb_compare_results
species_list=("s_cerevisiae" "c_elegans" "d_melanogaster" "a_thaliana" "m_musculus" "h_sapiens")
for species in "${species_list[@]}"
do
 compare_species_info "${species}"
 mkdir "${species}_CATS_rb_comparison" && cd "${species}_CATS_rb_comparison"
 for dir in $(find ../../ ../../../ref_transcriptomes/reference_trans_mapped -maxdepth 1 -type d -name *"${ref_transcriptome}"*)
 do
  find "${dir}" -type d -name "*CATS_rb_map" -exec ln -s {} . \;
 done
 for name in *_100_*
 do
  new_name="$(echo "${name}" | sed -E 's/.*100_([^/]*?)_CATS_rb_map$/\1/')"
  mv "${name}" "${new_name}"
 done
 for name in *_mult5_*
 do
  name_tmp="$(echo "${name}" | sed -E 's/.*(RSP|TRI)/\1/; s/_mult5_redundancy//; s/^[^_]*_[^_]*_[^_]*_//')"
  if [[ "$name_tmp" == *110000000 ]]
  then
   new_name="${name_tmp/%110000000/3}"
  elif [[ "$name_tmp" == *60000000 ]]
  then
   new_name="${name_tmp/%60000000/2}"
  elif [[ "$name_tmp" == *10000000 ]]
  then
   new_name="${name_tmp/%10000000/1}"
  fi
  mv "${name}" "${new_name}"
 done
 mv *CATS_rb_map ref

 echo "Running CATS-rb comparison on the complete dataset"
 CATS_rb_compare -p "${min_exon_id_prop}" -e 30 -i "${max_intron_len}" -l "${min_exon_set_len}" -L "${min_tr_set_len}" -m "${max_tr_set_len}" -d 600 -f "${num_longest_scaff}" -B "${num_genomic_bins}" -H 2000 -F "../../../ref_annotation/${ref_annotation}" -t 20 -D CATS_rb_comparison_sim1 "../../../ref_genomes/${ref_genome}" 1_4_0.005_RSP 1_4_0.01_RSP 1_4_0.02_RSP 5_10_0.005_RSP 5_10_0.01_RSP 5_10_0.02_RSP 11_20_0.005_RSP 11_20_0.01_RSP 11_20_0.02_RSP 21_30_0.005_RSP 21_30_0.01_RSP 21_30_0.02_RSP 31_40_0.005_RSP 31_40_0.01_RSP 31_40_0.02_RSP 41_50_0.005_RSP 41_50_0.01_RSP 41_50_0.02_RSP 51_100_0.005_RSP 51_100_0.01_RSP 51_100_0.02_RSP 1_4_0.005_TRI 1_4_0.01_TRI 1_4_0.02_TRI 5_10_0.005_TRI 5_10_0.01_TRI 5_10_0.02_TRI 11_20_0.005_TRI 11_20_0.01_TRI 11_20_0.02_TRI 21_30_0.005_TRI 21_30_0.01_TRI 21_30_0.02_TRI 31_40_0.005_TRI 31_40_0.01_TRI 31_40_0.02_TRI 41_50_0.005_TRI 41_50_0.01_TRI 41_50_0.02_TRI 51_100_0.005_TRI 51_100_0.01_TRI 51_100_0.02_TRI 1_4_0.005_IDB 1_4_0.01_IDB 1_4_0.02_IDB 5_10_0.005_IDB 5_10_0.01_IDB 5_10_0.02_IDB 11_20_0.005_IDB 11_20_0.01_IDB 11_20_0.02_IDB 21_30_0.005_IDB 21_30_0.01_IDB 21_30_0.02_IDB 31_40_0.005_IDB 31_40_0.01_IDB 31_40_0.02_IDB 41_50_0.005_IDB 41_50_0.01_IDB 41_50_0.02_IDB 51_100_0.005_IDB 51_100_0.01_IDB 51_100_0.02_IDB 1_4_0.005_SOA 1_4_0.01_SOA 1_4_0.02_SOA 5_10_0.005_SOA 5_10_0.01_SOA 5_10_0.02_SOA 11_20_0.005_SOA 11_20_0.01_SOA 11_20_0.02_SOA 21_30_0.005_SOA 21_30_0.01_SOA 21_30_0.02_SOA 31_40_0.005_SOA 31_40_0.01_SOA 31_40_0.02_SOA 41_50_0.005_SOA 41_50_0.01_SOA 41_50_0.02_SOA 51_100_0.005_SOA 51_100_0.01_SOA 51_100_0.02_SOA ref
 cd CATS_rb_comparison_sim1 && mv CATS_rb_main_comparison_results.tsv "${species}_CATS_rb_sim_rel_results1" && mv CATS_rb_annotation_based_analysis_results.tsv "${species}_CATS_rb_sim_annot_results" && cd ..

 echo "Running CATS-rb comparison on subset 1 (51-100x)"
 CATS_rb_compare -p "${min_exon_id_prop}" -e 30 -i "${max_intron_len}" -l "${min_exon_set_len}" -L "${min_tr_set_len}" -m "${max_tr_set_len}" -d 600 -f "${num_longest_scaff}" -B "${num_genomic_bins}" -H 2000 -t 20 -D CATS_rb_comparison_sim2 "../../../ref_genomes/${ref_genome}" 1_4_0.005_RSP 1_4_0.01_RSP 1_4_0.02_RSP 5_10_0.005_RSP 5_10_0.01_RSP 5_10_0.02_RSP 11_20_0.005_RSP 11_20_0.01_RSP 11_20_0.02_RSP 21_30_0.005_RSP 21_30_0.01_RSP 21_30_0.02_RSP 31_40_0.005_RSP 31_40_0.01_RSP 31_40_0.02_RSP 41_50_0.005_RSP 41_50_0.01_RSP 41_50_0.02_RSP 51_100_0.005_RSP 51_100_0.01_RSP 51_100_0.02_RSP 1_4_0.005_TRI 1_4_0.01_TRI 1_4_0.02_TRI 5_10_0.005_TRI 5_10_0.01_TRI 5_10_0.02_TRI 11_20_0.005_TRI 11_20_0.01_TRI 11_20_0.02_TRI 21_30_0.005_TRI 21_30_0.01_TRI 21_30_0.02_TRI 31_40_0.005_TRI 31_40_0.01_TRI 31_40_0.02_TRI 41_50_0.005_TRI 41_50_0.01_TRI 41_50_0.02_TRI 51_100_0.005_TRI 51_100_0.01_TRI 51_100_0.02_TRI 1_4_0.005_IDB 1_4_0.01_IDB 1_4_0.02_IDB 5_10_0.005_IDB 5_10_0.01_IDB 5_10_0.02_IDB 11_20_0.005_IDB 11_20_0.01_IDB 11_20_0.02_IDB 21_30_0.005_IDB 21_30_0.01_IDB 21_30_0.02_IDB 31_40_0.005_IDB 31_40_0.01_IDB 31_40_0.02_IDB 41_50_0.005_IDB 41_50_0.01_IDB 41_50_0.02_IDB 51_100_0.005_IDB 51_100_0.01_IDB 51_100_0.02_IDB 1_4_0.005_SOA 1_4_0.01_SOA 1_4_0.02_SOA 5_10_0.005_SOA 5_10_0.01_SOA 5_10_0.02_SOA 11_20_0.005_SOA 11_20_0.01_SOA 11_20_0.02_SOA 21_30_0.005_SOA 21_30_0.01_SOA 21_30_0.02_SOA 31_40_0.005_SOA 31_40_0.01_SOA 31_40_0.02_SOA 41_50_0.005_SOA 41_50_0.01_SOA 41_50_0.02_SOA 51_100_0.005_SOA 51_100_0.01_SOA 51_100_0.02_SOA
 cd CATS_rb_comparison_sim2 && mv CATS_rb_main_comparison_results.tsv "${species}_CATS_rb_sim_rel_results2" && cd ..

 echo "Running CATS-rb comparison on subset 2 (31-40x)"
 CATS_rb_compare -p "${min_exon_id_prop}" -e 30 -i "${max_intron_len}" -l "${min_exon_set_len}" -L "${min_tr_set_len}" -m "${max_tr_set_len}" -d 600 -f "${num_longest_scaff}" -B "${num_genomic_bins}" -H 2000 -t 20 -D CATS_rb_comparison_sim3 "../../../ref_genomes/${ref_genome}" 1_4_0.005_RSP 1_4_0.01_RSP 1_4_0.02_RSP 5_10_0.005_RSP 5_10_0.01_RSP 5_10_0.02_RSP 11_20_0.005_RSP 11_20_0.01_RSP 11_20_0.02_RSP 21_30_0.005_RSP 21_30_0.01_RSP 21_30_0.02_RSP 31_40_0.005_RSP 31_40_0.01_RSP 31_40_0.02_RSP 1_4_0.005_TRI 1_4_0.01_TRI 1_4_0.02_TRI 5_10_0.005_TRI 5_10_0.01_TRI 5_10_0.02_TRI 11_20_0.005_TRI 11_20_0.01_TRI 11_20_0.02_TRI 21_30_0.005_TRI 21_30_0.01_TRI 21_30_0.02_TRI 31_40_0.005_TRI 31_40_0.01_TRI 31_40_0.02_TRI 1_4_0.005_IDB 1_4_0.01_IDB 1_4_0.02_IDB 5_10_0.005_IDB 5_10_0.01_IDB 5_10_0.02_IDB 11_20_0.005_IDB 11_20_0.01_IDB 11_20_0.02_IDB 21_30_0.005_IDB 21_30_0.01_IDB 21_30_0.02_IDB 31_40_0.005_IDB 31_40_0.01_IDB 31_40_0.02_IDB 1_4_0.005_SOA 1_4_0.01_SOA 1_4_0.02_SOA 5_10_0.005_SOA 5_10_0.01_SOA 5_10_0.02_SOA 11_20_0.005_SOA 11_20_0.01_SOA 11_20_0.02_SOA 21_30_0.005_SOA 21_30_0.01_SOA 21_30_0.02_SOA 31_40_0.005_SOA 31_40_0.01_SOA 31_40_0.02_SOA
 cd CATS_rb_comparison_sim3 && mv CATS_rb_main_comparison_results.tsv "${species}_CATS_rb_sim_rel_results3" && cd ..

 echo "Running CATS-rb comparison on subset 3 (11-20x)"
 CATS_rb_compare -p "${min_exon_id_prop}" -e 30 -i "${max_intron_len}" -l "${min_exon_set_len}" -L "${min_tr_set_len}" -m "${max_tr_set_len}" -d 600 -f "${num_longest_scaff}" -B "${num_genomic_bins}" -H 2000 -t 20 -D CATS_rb_comparison_sim4 "../../../ref_genomes/${ref_genome}" 1_4_0.005_RSP 1_4_0.01_RSP 1_4_0.02_RSP 5_10_0.005_RSP 5_10_0.01_RSP 5_10_0.02_RSP 11_20_0.005_RSP 11_20_0.01_RSP 11_20_0.02_RSP 1_4_0.005_TRI 1_4_0.01_TRI 1_4_0.02_TRI 5_10_0.005_TRI 5_10_0.01_TRI 5_10_0.02_TRI 11_20_0.005_TRI 11_20_0.01_TRI 11_20_0.02_TRI 1_4_0.005_IDB 1_4_0.01_IDB 1_4_0.02_IDB 5_10_0.005_IDB 5_10_0.01_IDB 5_10_0.02_IDB 11_20_0.005_IDB 11_20_0.01_IDB 11_20_0.02_IDB 1_4_0.005_SOA 1_4_0.01_SOA 1_4_0.02_SOA 5_10_0.005_SOA 5_10_0.01_SOA 5_10_0.02_SOA 11_20_0.005_SOA 11_20_0.01_SOA 11_20_0.02_SOA
 cd CATS_rb_comparison_sim4 && mv CATS_rb_main_comparison_results.tsv "${species}_CATS_rb_sim_rel_results4" && cd ..

 echo "Running CATS-rb comparison on subset 4 (5-10x)"
 CATS_rb_compare -p "${min_exon_id_prop}" -e 30 -i "${max_intron_len}" -l "${min_exon_set_len}" -L "${min_tr_set_len}" -m "${max_tr_set_len}" -d 600 -f "${num_longest_scaff}" -B "${num_genomic_bins}" -H 2000 -t 20 -D CATS_rb_comparison_sim5 "../../../ref_genomes/${ref_genome}" 1_4_0.005_RSP 1_4_0.01_RSP 1_4_0.02_RSP 5_10_0.005_RSP 5_10_0.01_RSP 5_10_0.02_RSP 1_4_0.005_TRI 1_4_0.01_TRI 1_4_0.02_TRI 5_10_0.005_TRI 5_10_0.01_TRI 5_10_0.02_TRI 1_4_0.005_IDB 1_4_0.01_IDB 1_4_0.02_IDB 5_10_0.005_IDB 5_10_0.01_IDB 5_10_0.02_IDB 1_4_0.005_SOA 1_4_0.01_SOA 1_4_0.02_SOA 5_10_0.005_SOA 5_10_0.01_SOA 5_10_0.02_SOA
 cd CATS_rb_comparison_sim5 && mv CATS_rb_main_comparison_results.tsv "${species}_CATS_rb_sim_rel_results5" && cd ..

 echo "Running CATS-rb comparison on subset 5 (1-4x)"
 CATS_rb_compare -p "${min_exon_id_prop}" -e 30 -i "${max_intron_len}" -l "${min_exon_set_len}" -L "${min_tr_set_len}" -m "${max_tr_set_len}" -d 600 -f "${num_longest_scaff}" -B "${num_genomic_bins}" -H 2000 -t 20 -D CATS_rb_comparison_sim6 "../../../ref_genomes/${ref_genome}" 1_4_0.005_RSP 1_4_0.01_RSP 1_4_0.02_RSP 1_4_0.005_TRI 1_4_0.01_TRI 1_4_0.02_TRI 1_4_0.005_IDB 1_4_0.01_IDB 1_4_0.02_IDB 1_4_0.005_SOA 1_4_0.01_SOA 1_4_0.02_SOA
 cd CATS_rb_comparison_sim6 && mv CATS_rb_main_comparison_results.tsv "${species}_CATS_rb_sim_rel_results6" && cd ..

 echo "Running CATS-rb comparison on assemblies undergoing multiplicative mutation simulation"
 CATS_rb_compare -p "${min_exon_id_prop}" -e 30 -i "${max_intron_len}" -l "${min_exon_set_len}" -L "${min_tr_set_len}" -m "${max_tr_set_len}" -d 600 -f "${num_longest_scaff}" -B "${num_genomic_bins}" -H 2000 -t 20 -D CATS_rb_comparison_sim_mut "../../../ref_genomes/${ref_genome}" 21_30_0.005_RSP 21_30_0.005_TRI RSP_0.1_1 RSP_0.1_2 RSP_0.1_3 RSP_0.2_1 RSP_0.2_2 RSP_0.2_3 RSP_0.3_1 RSP_0.3_2 RSP_0.3_3 RSP_0.4_1 RSP_0.4_2 RSP_0.4_3 RSP_0.5_1 RSP_0.5_2 RSP_0.5_3 RSP_0.6_1 RSP_0.6_2 RSP_0.6_3 TRI_0.1_1 TRI_0.1_2 TRI_0.1_3 TRI_0.2_1 TRI_0.2_2 TRI_0.2_3 TRI_0.3_1 TRI_0.3_2 TRI_0.3_3 TRI_0.4_1 TRI_0.4_2 TRI_0.4_3 TRI_0.5_1 TRI_0.5_2 TRI_0.5_3 TRI_0.6_1 TRI_0.6_2 TRI_0.6_3
 cd CATS_rb_comparison_sim_mut && mv CATS_rb_main_comparison_results.tsv "${species}_CATS_rb_sim_mut_rel_results" && cd ..
 cd ..
done
cd ../../

#Mapping and analysing reference transcriptomes undergoing chimerism simulation
cd ref_transcriptomes/chimerism
chim_ref_transcriptomes="$(ls -tr *chimerism*)"
for chim_ref_transcriptome in ${chim_ref_transcriptomes}
do
 map_species_info "${chim_ref_transcriptome}"
 compare_species_info "${chim_ref_transcriptome}"
 cat "${chim_ref_transcriptome}" | grep ">" | tr -d '>' | awk '{print $1}' > "${chim_ref_transcriptome}_names"
 echo "Mapping and analysing ${chim_ref_transcriptome}"
 CATS_rb_map -p "${species_preset}" -t 20 "../../ref_genomes/ref_genome_indices/${ref_genome_index}" "${chim_ref_transcriptome}"
 CATS_rb_compare -p "${min_exon_id_prop}" -e 30 -i "${max_intron_len}" -d 600 -f "${num_longest_scaff}" -B "${num_genomic_bins}" -t 20 -D "${chim_ref_transcriptome}_CATS_rb_comparison" "../../ref_genomes/${ref_genome}" "${chim_ref_transcriptome}_CATS_rb_map"
 cd "${chim_ref_transcriptome}_CATS_rb_comparison" && mv str_inconsistent_transcripts.tsv "../${chim_ref_transcriptome}_str_inconsistent_transcripts.tsv" && cd ..
done
cd ../../

#Downloading public RNA-seq data from SRA
sra_ids=(SRR32108057 SRR30685385 SRR27822254 SRR26147123 SRR26001199 SRR24356101 SRR10815431
SRR32732688 SRR32058534 SRR31107386 SRR12868377 SRR17736005 SRR5224028 SRR6266308
SRR31530981 SRR29606991 SRR24876155 SRR23933588 SRR23870057 SRR20326862 SRR14560308
SRR33339778 SRR30855151 SRR14160649 SRR24948904 SRR24636232 SRR13159213 SRR8422221
SRR32903297 SRR31588239 SRR27369179 SRR28908270 SRR18855849 SRR24134730 SRR13981556
SRR32139880 SRR31785722 SRR25646397 SRR25732618 SRR22548603 SRR8357441 SRR7741229)

mkdir public_data && cd public_data
for sra_id in "${sra_ids[@]}"
do
 case "${sra_id}" in
  SRR32108057|SRR30685385|SRR27822254|SRR26147123|SRR26001199|SRR24356101|SRR10815431) species=s_cerevisiae ;;
  SRR32732688|SRR32058534|SRR31107386|SRR12868377|SRR17736005|SRR5224028|SRR6266308)   species=c_elegans ;;
  SRR31530981|SRR29606991|SRR24876155|SRR23933588|SRR23870057|SRR20326862|SRR14560308) species=d_melanogaster ;;
  SRR33339778|SRR30855151|SRR14160649|SRR24948904|SRR24636232|SRR13159213|SRR8422221)  species=a_thaliana ;;
  SRR32903297|SRR31588239|SRR27369179|SRR28908270|SRR18855849|SRR24134730|SRR13981556) species=m_musculus ;;
  SRR32139880|SRR31785722|SRR25646397|SRR25732618|SRR22548603|SRR8357441|SRR7741229)   species=h_sapiens ;;
 esac
 echo "Downloading ${sra_id}"
 prefetch -v "${sra_id}"
 fastq-dump --split-files -O "pub_${species}_${sra_id}" "${sra_id}"
 cd "pub_${species}_${sra_id}"
 mv *_1.fastq reads1.fq
 mv *_2.fastq reads2.fq
 cd ..
done

#Trimming adapters and low-quality bases in public libraries
for dir in pub_*
do
 cd "${dir}"
 trim_galore --paired --phred33 -q 5 --stringency 5 -j 4 reads1.fq reads2.fq
 cd ..
done

#Assembling public libraries with rnaSPAdes, Trinity, IDBA-tran, and SOAPdenovo-Trans
for dir in pub_*
do
 cd "${dir}"

 outdir="${dir}_rnaspades_dir"
 echo "Assembling ${dir} with rnaSPAdes"
 rnaspades.py -1 reads1_val_1.fq -2 reads2_val_2.fq -t 20 -m 300 -o "${outdir}"
 mv "${outdir}/transcripts.fasta" "${dir}_RSP"

 outdir="${dir}_trinity_dir"
 echo "Assembling ${dir} with Trinity"
 Trinity --left reads1_val_1.fq --right reads2_val_2.fq --seqType fq --CPU 20 --max_memory 300G --output "${outdir}"
 mv "${outdir}.Trinity.fasta" "${dir}_TRI"

 outdir="${dir}_idba_dir"
 fq2fa --merge --filter reads1_val_1.fq reads2_val_2.fq reads_together.fq
 echo "Assembling ${dir} with IDBA-tran"
 idba_tran -r reads_together.fq --num_threads 20 -o "${outdir}"
 mv "${outdir}/contig.fa" "${dir}_IDB"

 READ_LEN="$(head -n 4000000 reads1_val_1.fq | awk 'NR%4 == 2 {print length}' | awk '{ total += $1 } END { print total / NR }')"
 READ_LEN="$(printf "%.0f" "${READ_LEN}")"
 echo "Assembling ${dir} with SOAPdenovo-Trans"
 SOAPdenovo-Trans-31mer all -s ../../soap_configfile_public -p 20 -o "${dir}_soap"
 mv "${dir}_soap.scafSeq" "${dir}_SOA"
 cd ..
done

#Running CATS-rf, TransRate, and RSEM-EVAL on public assemblies
for dir in pub_*
do
 cd "${dir}"
 for transcriptome in *_RSP *_TRI *_IDB *_SOA
 do

  echo "Running CATS-rf on ${transcriptome}"
  CATS_rf -t 20 -M 16384M -o "${transcriptome}_CATS_rf" "${transcriptome}" reads1_val_1.fq reads2_val_2.fq

  mkdir "${transcriptome}_rsem_eval" && cd "${transcriptome}_rsem_eval"
  echo "Running RSEM-EVAL on ${transcriptome}"
  rsem-eval-estimate-transcript-length-distribution "../${transcriptome}" parameter_file
  kallisto index -t 20 -i "${transcriptome}_kallisto_idx" "../${transcriptome}"
  kallisto quant -i "${transcriptome}_kallisto_idx" -t 20 -o "${transcriptome}_kallisto_dir" ../reads1_val_1.fq ../reads2_val_2.fq 2> "${transcriptome}_kallisto_log"
  cat "${transcriptome}_kallisto_log" | grep "fragment length:" | cut -d ':' -f2- | xargs > frag_length_file
  FRAG_LENGTH="$(< frag_length_file)"
  rsem-eval-calculate-score --seed 100 -p 20 --paired-end --transcript-length-parameters parameter_file ../reads1_val_1.fq ../reads2_val_2.fq "../${transcriptome}" "${transcriptome}_rsem" "${FRAG_LENGTH}"
  cut -f 2 "${transcriptome}_rsem.score" | head -n 1 > "${transcriptome}_rsem_eval_assembly_score"
  cut -f 1,9 "${transcriptome}_rsem.score.isoforms.results" > "${transcriptome}_rsem_eval_transcript_scores"
  cd ..

  echo "Running TransRate on ${transcriptome}"
  transrate --assembly="${transcriptome}" --left=reads1_val_1.fq --right=reads2_val_2.fq --threads=20 --output="${transcriptome}_transrate"
  cut -d "," -f 37 "${transcriptome}_transrate/assemblies.csv" | tail -n +2 > "${transcriptome}_transrate/${transcriptome}_transrate_assembly_score"
  cut -d "," -f 1,16,15,17,18,9 "${transcriptome}"_transrate/pub*/contigs.csv > "${transcriptome}_transrate/${transcriptome}_transrate_transcript_scores"
 done
 cd ..
done

#Calculating F-scores of public assemblies with CRB-BLAST
for dir in pub_*
do
 cd "${dir}"
 map_ref_transcriptome_info "${dir}"
 for transcriptome in *_RSP *_TRI *_IDB *_SOA
 do
  mkdir "${transcriptome}_crb" && cd "${transcriptome}_crb"
  echo "Running BLAT on ${transcriptome}"
  pblat "../../../ref_transcriptomes/${ref_transcriptome}" "../${transcriptome}" -t=dna -q=dna -maxIntron=0 -tileSize=8 -threads=20 -out=psl blat_output_tmp.psl
  tail -n +6 blat_output_tmp.psl > blat_output.psl
  filter_by_blat.R 20 blat_output.psl 0.5 "../${transcriptome}" "${transcriptome}"
  echo "Running CRB-BLAST on ${transcriptome}"
  crb-blast -q "filtered_0.5_${transcriptome}" -t "../../../ref_transcriptomes/${ref_transcriptome}" -h 20 -o "${transcriptome}_crb"
  get_f_scores_from_crb_table.R 20 "${transcriptome}_crb" "filtered_0.5_${transcriptome}" "filtered_0.5_${transcriptome}"
  cd ..
 done
 cd ..
done

#Mapping public transcriptome assemblies to the respective reference genomes with CATS_rb_map
for dir in pub_*
do
 cd "${dir}"
 map_species_info "${dir}"
 for transcriptome in *_RSP *_TRI *_IDB *_SOA
 do
  echo "Mapping ${transcriptome} to the reference genome"
  CATS_rb_map -p "${species_preset}" -t 20 "../../ref_genomes/ref_genome_indices/${ref_genome_index}" "${transcriptome}"
 done
 cd ..
done

#Comparing public transcriptome assemblies with CATS_rb_compare
mkdir CATS_rb_compare_results && cd CATS_rb_compare_results
for species in "${species_list[@]}"
do
 compare_species_info "${species}"
 mkdir "${species}_CATS_rb_comparison" && cd "${species}_CATS_rb_comparison"
 for dir in $(find ../../ ../../../ref_transcriptomes/reference_trans_mapped -maxdepth 1 -type d -name *"${species}"*)
 do
  find "${dir}" -type d -name "*CATS_rb_map" -exec ln -s {} . \;
 done
 for name in pub*_SRR*
 do
  new_name="$(echo "${name}" | sed -E 's/.*(SRR[0-9]+)_(SOA|TRI|RSP|IDB)_CATS_rb_map/\1_\2/')"
  mv "${name}" "${new_name}"
 done
 mv *CATS_rb_map ref
 sra_ids="$(ls | grep SRR | sed 's/_.*//')"
 for sra_id in "${sra_ids[@]}"
 do
  CATS_rb_compare -p "${min_exon_id_prop}" -e 30 -i "${max_intron_len}" -l "${min_exon_set_len}" -L "${min_tr_set_len}" -m "${max_tr_set_len}" -d 600 -f "${num_longest_scaff}" -B "${num_genomic_bins}" -H 2000 -F "../../../ref_annotation/${ref_annotation}" -t 20 -D CATS_rb_comparison_pub "../../../ref_genomes/${ref_genome}" "${sra_id}_RSP" "${sra_id}_TRI" "${sra_id}_IDB" "${sra_id}_SOA" ref
  cd CATS_rb_comparison_pub && mv CATS_rb_main_comparison_results.tsv "${species}_${sra_id}_CATS_rb_pub_rel_results" && mv CATS_rb_annotation_based_analysis_results.tsv "${species}_${sra_id}_CATS_rb_pub_annot_results" && cd ..
 done
 cd ..
done
cd ../../

#Preparing for R analysis: CATS-rf simulated and public assemblies (Figure 2)
mkdir figure_2_data && cd figure_2_data

#Linking CATS-rf transcript scores of simulated assemblies
find ../simulated_data -type f -name "sim_*_CATS_rf_transcript_scores.tsv" -exec ln -s {} . \;

#Linking RSEM-EVAL transcript scores of simulated assemblies
find ../simulated_data -type f -name "sim_*_rsem_eval_transcript_scores" -exec ln -s {} . \;

#Linking TransRate transcript scores of simulated assemblies
find ../simulated_data -type f -name "sim_*_transrate_transcript_scores" -exec ln -s {} . \;

#Linking RSEM-EVAL assembly scores of simulated assemblies
find ../simulated_data -type f -name "sim_*_rsem_eval_assembly_score" -exec ln -s {} . \;

#Linking TransRate assembly scores of simulated assemblies
find ../simulated_data -type f -name "sim_*_transrate_assembly_score" -exec ln -s {} . \;

#Linking transcript F-scores of simulated assemblies
find ../simulated_data -type f -name "sim_*_f_scores" -exec ln -s {} . \;

#Linking CATS-rf transcript scores of public assemblies
find ../public_data -type f -name "pub*_*_CATS_rf_transcript_scores.tsv" -exec ln -s {} . \;

#Linking RSEM-EVAL transcript scores of public assemblies
find ../public_data -type f -name "pub_*_rsem_eval_transcript_scores" -exec ln -s {} . \;

#Linking TransRate transcript scores of public assemblies
find ../public_data -type f -name "pub_*_transrate_transcript_scores" -exec ln -s {} . \;

#Linking transcript F-scores of public assemblies
find ../public_data -type f -name "pub_*_f_scores" -exec ln -s {} . \;

#Combining the linked CATS-rf results into a single table
merge_cats_rf_results.R

#Generating Figure 2 elements
generate_fig2_elements.R
cd ..

#Preparing for R analysis: CATS-rf mutated simulated assemblies (Figure 3)
mkdir figure_3_data && cd figure_3_data

#Linking CATS-rf transcript scores of mutated simulated assemblies. Merging transcript scores into a single table
find ../simulated_data -type f -name "sim_*_21_30_0.005_{RSP,TRI}_CATS_rf_mut_*_transcript_scores.tsv" -exec ln -s {} . \;
echo -e "transcript\tcats_rf_coverage_component\tcats_rf_accuracy_component\tcats_rf_local_fidelity_component\tcats_rf_integrity_component\tcats_rf_transcript_score\tassembly" > simulated_mutated_assemblies_CATS_rf_transcript_scores
for file in sim_*_21_30_0.005_{RSP,TRI}_CATS_rf_mut_*_transcript_scores.tsv
do
 tail -n +2 "${file}" | awk -v fname="${file}" -F'\t' 'BEGIN{OFS="\t"} {print $0, fname}' >> simulated_mutated_assemblies_CATS_rf_transcript_scores
done

find ../simulated_data -type f -name "sim_*_21_30_0.005_{RSP,TRI}_CATS_rf_transcript_scores.tsv" -exec ln -s {} . \;
for file in sim_*_21_30_0.005_{RSP,TRI}_CATS_rf_transcript_scores.tsv
do
 tail -n +2 "${file}" | awk -v fname="${file}" -F'\t' 'BEGIN{OFS="\t"} {print $0, fname}' >> simulated_native_assemblies_CATS_rf_transcript_scores
done

cat simulated_mutated_assemblies_CATS_rf_transcript_scores simulated_native_assemblies_CATS_rf_transcript_scores | cut -f 7,1,2,3,4,5,6 > mutation_analysis_CATS_rf_transcript_scores_for_figure3.tsv

#Generating Figure 3 elements
generate_fig3_elements.R
cd ..

#Preparing for R analysis: CATS-rb simulated and public assemblies (Figure 5)
mkdir figure_5_data && cd figure_5_data

#Linking CATS-rb scores of simulated assemblies
find ../simulated_data -type f -name "*_CATS_rb_sim_rel_results*" -exec ln -s {} . \;
find ../simulated_data -type f -name "*_CATS_rb_sim_annot_results" -exec ln -s {} . \;

#Linking CATS-rb scores of mutated simulated assemblies
find ../simulated_data -type f -name "*_CATS_rb_sim_mut_rel_results" -exec ln -s {} . \;

#Linking CATS-rb scores of public assemblies
find ../public_data -type f -name "*_CATS_rb_pub_rel_results" -exec ln -s {} . \;
find ../public_data -type f -name "*_CATS_rb_pub_annot_results" -exec ln -s {} . \;

#Extracting exon and transcript scores
for file in *rel_results*
do
 sed -n '1p;18p;31p' "${file}" | cut -f2- > "${file}.tmp" && mv "${file}.tmp" "${file}"
done

for file in *annot_results*
do
 sed -n '1p;6p;11p' "${file}" | cut -f2- > "${file}.tmp" && mv "${file}.tmp" "${file}"
done

#Renaming public assemblies
HEADER="REF\tRSP\tTRI\tIDB\tSOA"
for file in *_CATS_rb_pub_*
do
 sed -i "1s/.*/${HEADER}/" "${file}"
done

#Combining the linked CATS-rb results into a single table
merge_cats_rb_results.R

#Linking CATS-rf assembly score and transcript F-score mean table (from Figure 2)
ln -s ../figure_2_data/simulated_cats_rf_f_results_for_figure5.tsv

#Generating Figure 5 elements
generate_fig5_elements.R
cd ..

#Preparing for R analysis: CATS-rb reference chimeric transcripts (Extended data figure 7)
mkdir extended_data_figure_7_data && cd extended_data_figure_7_data
find ../ref_transcriptomes/chimerism -type f -name "*_chimerism_*names" -exec ln -s {} . \;
find ../ref_transcriptomes/chimerism -type f -name "*str_inconsistent_transcripts.tsv" -exec ln -s {} . \;

#Merging transcript names into a single table
echo -e "transcript\tsource_file" > chimeric_ref_transcriptome_names_for_figure_ext7
for file in *_chimerism_*names
do
 cat "${file}" | awk -v fname="${file}" -F'\t' 'BEGIN{OFS="\t"} {print $0, fname}' >> chimeric_ref_transcriptome_names_for_figure_ext7
done

#Merging structurally inconsistent transcript lists into a single table
echo -e "transcript\tsource_file" > str_inconsistent_transcripts_for_figure_ext7
for file in *str_inconsistent_transcripts.tsv
do
 cut -f 1 "${file}" | tail -n +2 | awk -v fname="${file}" -F'\t' 'BEGIN{OFS="\t"} {print $0, fname}' >> str_inconsistent_transcripts_for_figure_ext7
done

#Generating Extended data figure 7
generate_fig_ext7.R
cd ..

exit 0
