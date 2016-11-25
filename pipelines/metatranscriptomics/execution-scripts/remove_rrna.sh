#!/bin/bash
#PBS -l vmem=40g
#PBS -l nodes=1:ppn=6
#PBS -l walltime=48:00:00
#PBS -d /hpf/largeprojects/ccmbio/nreinhardt/mugqic_pipelines/pipelines/metatranscriptomics/execution-scripts

module load mugqic-pipelines/2.2.0
module load infernal
module load perl

rfam_path=/hpf/largeprojects/ccmbio/nreinhardt/metatranscriptomics-test/reference-files/Rfam.cm

input_dir=remove_duplicates
output_dir=remove_rrna

mkdir $output_dir

cmscan -o $output_dir/cow1_rRNA.log --tblout $output_dir/cow1_rRNA.infernalout --noali --notextw --rfam -E 0.001 $rfam_path $input_dir/cow1_qual_all_unique.fasta
cmscan -o $output_dir/cow2_rRNA.log --tblout $output_dir/cow2_rRNA.infernalout --noali --notextw --rfam -E 0.001 $rfam_path $input_dir/cow2_qual_all_unique.fasta


# cow1
python /hpf/largeprojects/ccmbio/nreinhardt/mugqic_pipelines/pipelines/metatranscriptomics/scripts/main_get_sequence_length.py \
    --fasta $input_dir/cow1_qual_all_unique.fasta \
    --id-file $input_dir/cow1_qual_all_unique_IDs.txt \
    --output $output_dir/cow1_IDs_length.txt
# cow2
python /hpf/largeprojects/ccmbio/nreinhardt/mugqic_pipelines/pipelines/metatranscriptomics/scripts/main_get_sequence_length.py \
    --fasta $input_dir/cow2_qual_all_unique.fasta \
    --id-file $input_dir/cow2_qual_all_unique_IDs.txt \
    --output $output_dir/cow2_IDs_length.txt
#python /hpf/largeprojects/ccmbio/nreinhardt/mugqic_pipelines/pipelines/metatranscriptomics/scripts/main_get_sequence_length.py
#perl main_get_sequence_length.pl cow rRNA

perl /hpf/largeprojects/ccmbio/nreinhardt/mugqic_pipelines/pipelines/metatranscriptomics/scripts/main_get_infernal_fromfile_1tophit.pl cow pairs 1 0.001 90
#perl main_get_infernal_fromfile_1tophit.pl cow pairs 1 0.001 90

# cow1
python /hpf/largeprojects/ccmbio/nreinhardt/mugqic_pipelines/pipelines/metatranscriptomics/scripts/split_reads_by_id.py \
    --fastq $input_dir/cow1_qual_all_unique.fastq \
    --id-file $output_dir/cow1_rRNA_infernal_IDs.txt \
    --included $output_dir/cow1_qual_unique_rRNA.fastq \
    --excluded $output_dir/cow1_qual_unique_n_rRNA.fastq
# cow2
python /hpf/largeprojects/ccmbio/nreinhardt/mugqic_pipelines/pipelines/metatranscriptomics/scripts/split_reads_by_id.py \
    --fastq $input_dir/cow2_qual_all_unique.fastq \
    --id-file $output_dir/cow2_rRNA_infernal_IDs.txt \
    --included $output_dir/cow2_qual_unique_rRNA.fastq \
    --excluded $output_dir/cow2_qual_unique_n_rRNA.fastq
#perl /hpf/largeprojects/ccmbio/nreinhardt/mugqic_pipelines/pipelines/metatranscriptomics/scripts/main_select_reads_fromfile.pl
#perl main_select_reads_fromfile.pl cow rRNA infernal pairs
