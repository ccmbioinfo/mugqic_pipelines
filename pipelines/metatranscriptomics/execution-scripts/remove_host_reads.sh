#!/bin/bash

module load mugqic-pipelines/2.2.0
#module load mugqic/bwa/0.7.12
module load perl

bwa_dir=../reference-files/bwa-0.7.5a
input_dir=remove_rrna
output_dir=remove_host_reads

mkdir $output_dir

COW_CDS=/hpf/largeprojects/ccmbio/nreinhardt/metatranscriptomics-test/reference-files/cow_cds.fa


# Keep only reads in both cow1 and cow2
#python /hpf/largeprojects/ccmbio/nreinhardt/mugqic_pipelines/pipelines/metatranscriptomics/scripts/fastq_intersection.py \
#    --fastq1 $input_dir/cow1_qual_unique_n_rRNA.fastq \
#    --fastq2 $input_dir/cow2_qual_unique_n_rRNA.fastq \
#    --output1 $output_dir/cow1_matching_ids.fastq \
#    --output2 $output_dir/cow2_matching_ids.fastq

$bwa_dir/bwa aln -t 4 $COW_CDS $input_dir/cow1_qual_unique_n_rRNA.fastq > $output_dir/cow1_host.sai
$bwa_dir/bwa aln -t 4 $COW_CDS $input_dir/cow2_qual_unique_n_rRNA.fastq > $output_dir/cow2_host.sai

#$bwa_dir/bwa sampe $COW_CDS $output_dir/cow1_host.sai $output_dir/cow2_host.sai $output_dir/cow1_matching_ids.fastq $output_dir/cow2_matching_ids.fastq > $output_dir/cow_host.sam
$bwa_dir/bwa sampe $COW_CDS \
    $output_dir/cow1_host.sai $output_dir/cow2_host.sai \
    $input_dir/cow1_qual_unique_n_rRNA.fastq $input_dir/cow2_qual_unique_n_rRNA.fastq \
    > $output_dir/cow_host.sam
#$bwa_dir/bwa sampe cow_cds.fa cow1_host.sai cow2_host.sai cow1_qual_unique_n_rRNA.fastq cow2_qual_unique_n_rRNA.fastq > cow_host.sam

# extract unmapped reads
samtools view -bS $output_dir/cow_host.sam | samtools sort -n -o $output_dir/cow_host.bam
samtools view -F 4 $output_dir/cow_host.bam > $output_dir/cow_host.bwaout

perl /hpf/largeprojects/ccmbio/nreinhardt/mugqic_pipelines/pipelines/metatranscriptomics/scripts/main_read_samout.pl \
    $output_dir/cow_host.bwaout \
    $output_dir/cow_host_bwa_IDs.txt \
    $output_dir/cow_host_bwa_pairs.txt \
    $output_dir/cow_host_bwa_hitsID.txt
#perl /hpf/largeprojects/ccmbio/nreinhardt/mugqic_pipelines/pipelines/metatranscriptomics/scripts/main_read_samout.pl cow host bwa pairs

# cow1
/hpf/largeprojects/ccmbio/nreinhardt/mugqic_pipelines/pipelines/metatranscriptomics/scripts/split_reads_by_id.py \
    --fastq $input_dir/cow1_qual_unique_n_rRNA.fastq \
    --id-file $output_dir/cow_host_bwa_IDs.txt \
    --included $output_dir/cow1_qual_unique_n_rRNA_host.fastq \
    --excluded $output_dir/cow1_qual_unique_n_rRNA_n_host.fastq
# cow2
/hpf/largeprojects/ccmbio/nreinhardt/mugqic_pipelines/pipelines/metatranscriptomics/scripts/split_reads_by_id.py \
    --fastq $input_dir/cow2_qual_unique_n_rRNA.fastq \
    --id-file $output_dir/cow_host_bwa_IDs.txt \
    --included $output_dir/cow2_qual_unique_n_rRNA_host.fastq \
    --excluded $output_dir/cow2_qual_unique_n_rRNA_n_host.fastq
#perl /hpf/largeprojects/ccmbio/nreinhardt/mugqic_pipelines/pipelines/metatranscriptomics/scripts/main_select_reads_fromfile.pl cow host bwa pairs
