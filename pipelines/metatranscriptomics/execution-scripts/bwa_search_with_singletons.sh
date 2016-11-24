#!/usr/bin/env bash

module load mugqic-pipelines/2.2.0
module load perl
module load mugqic/bwa/0.7.12
module load mugqic/samtools/1.3
module load blastp/2.2.27

input_dir=map_reads
output_dir=bwa_search_with_singletons
blastdb=../reference-files

mkdir $output_dir


bwa aln -t 4 $blastdb/microbial_all_cds.fasta $input_dir/cow1_singletons.fastq > $output_dir/cow1_singletons.sai
bwa aln -t 4 $blastdb/microbial_all_cds.fasta $input_dir/cow2_singletons.fastq > $output_dir/cow2_singletons.sai

python /hpf/largeprojects/ccmbio/nreinhardt/mugqic_pipelines/pipelines/metatranscriptomics/scripts/fastq_intersection.py \
    --fastq1 $input_dir/cow1_singletons.fastq \
    --fastq2 $input_dir/cow2_singletons.fastq \
    --output1 $input_dir/cow1_singletons.fastq \
    --output2 $input_dir/cow2_singletons.fastq

bwa sampe $blastdb/microbial_all_cds.fasta $output_dir/cow1_singletons.sai  $output_dir/cow2_singletons.sai $input_dir/cow1_singletons.fastq $input_dir/cow2_singletons.fastq > $output_dir/cow_singletons.sam
samtools view -bS $output_dir/cow_singletons.sam | samtools sort -n -o $output_dir/cow_singletons.bam
samtools view -F 4 $output_dir/cow_singletons.bam > $output_dir/cow_singletons_micro_cds.bwaout

perl /hpf/largeprojects/ccmbio/nreinhardt/mugqic_pipelines/pipelines/metatranscriptomics/scripts/main_read_samout.pl \
    $output_dir/cow_singletons_micro_cds.bwaout \
    $output_dir/cow_singletons_micro_cds_bwa_IDs.txt \
    $output_dir/cow_singletons_micro_cds_bwa_pairs.txt \
    $output_dir/cow_singletons_micro_cds_bwa_hitsID.txt
#perl main_read_samout.pl cow microgenes bwa singletons micro_cds


# cow1
python /hpf/largeprojects/ccmbio/nreinhardt/mugqic_pipelines/pipelines/metatranscriptomics/scripts/split_reads_by_id.py \
    --fastq $input_dir/cow1_singletons.fastq \
    --id-file $output_dir/cow_singletons_micro_cds_bwa_IDs.txt \
    --included $output_dir/cow1_singletons_micro_cds.fasta \
    --excluded $output_dir/cow1_singletons_n_micro_cds.fasta
# cow2
python /hpf/largeprojects/ccmbio/nreinhardt/mugqic_pipelines/pipelines/metatranscriptomics/scripts/split_reads_by_id.py \
    --fastq $input_dir/cow2_singletons.fastq \
    --id-file $output_dir/cow_singletons_micro_cds_bwa_IDs.txt \
    --included $output_dir/cow2_singletons_micro_cds.fasta \
    --excluded $output_dir/cow2_singletons_n_micro_cds.fasta
#perl main_select_reads_fromfile.pl cow microgenes bwa singletons micro_cds

