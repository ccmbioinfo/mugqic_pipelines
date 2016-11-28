#!/bin/bash

module load mugqic-pipelines/2.2.0
module load mugqic/bwa/0.7.12
module load perl

mkdir remove_host_reads

# Take only those reads that have matching ids in cow1 and cow2
#intersection-of-non-rnas.py
COW_CDS=/hpf/largeprojects/ccmbio/nreinhardt/metatranscriptomics-test/reference-files/cow_cds.fa

bwa aln -t 4 $COW_CDS remove_rrna/cow1_qual_unique_n_rRNA.fastq > remove_host_reads/cow1_host.sai
bwa aln -t 4 $COW_CDS remove_rrna/cow2_qual_unique_n_rRNA.fastq > remove_host_reads/cow2_host.sai
#bwa aln -t 4 cow_cds.fa cow1_qual_unique_n_rRNA.fastq > cow1_host.sai
#bwa aln -t 4 cow_cds.fa cow2_qual_unique_n_rRNA.fastq > cow2_host.sai

# Sort the fastqs by id
#/hpf/largeprojects/ccmbio/nreinhardt/mugqic_pipelines/pipelines/metatranscriptomics/scripts/sort_by_id.sh \
#    remove_rrna/cow1_qual_unique_n_rRNA.fastq
#/hpf/largeprojects/ccmbio/nreinhardt/mugqic_pipelines/pipelines/metatranscriptomics/scripts/sort_by_id.sh \
#    remove_rrna/cow2_qual_unique_n_rRNA.fastq

# Keep only reads in both cow1 and cow2
python /hpf/largeprojects/ccmbio/nreinhardt/mugqic_pipelines/pipelines/metatranscriptomics/scripts/fastq_intersection.py \
    --fastq1 remove_rrna/cow1_qual_unique_n_rRNA.fastq \
    --fastq2 remove_rrna/cow2_qual_unique_n_rRNA.fastq \
    --output1 remove_host_reads/cow1_matching_ids.fastq \
    --output2 remove_host_reads/cow2_matching_ids.fastq

bwa sampe $COW_CDS remove_host_reads/cow1_host.sai remove_host_reads/cow2_host.sai remove_host_reads/cow1_matching_ids.fastq remove_host_reads/cow2_matching_ids.fastq > remove_host_reads/cow_host.sam
#bwa sampe $COW_CDS remove_host_reads/cow1_host.sai remove_host_reads/cow2_host.sai remove_rrna/cow1_qual_unique_n_rRNA.fastq remove_rrna/cow2_qual_unique_n_rRNA.fastq > remove_host_reads/cow_host.sam
#bwa sampe cow_cds.fa cow1_host.sai cow2_host.sai cow1_qual_unique_n_rRNA.fastq cow2_qual_unique_n_rRNA.fastq > cow_host.sam

# extract unmapped reads
samtools view -bS remove_host_reads/cow_host.sam | samtools sort -n -o remove_host_reads/cow_host.bam
samtools view remove_host_reads/cow_host.bam > remove_host_reads/cow_host.bwaout
#samtools view -F 4 remove_host_reads/cow_host.bam > remove_host_reads/cow_host.bwaout

perl /hpf/largeprojects/ccmbio/nreinhardt/mugqic_pipelines/pipelines/metatranscriptomics/scripts/main_read_samout.pl \
    remove_host_reads/cow_host.bwaout \
    remove_host_reads/cow_host_bwa_IDs.txt \
    remove_host_reads/cow_host_bwa_pairs.txt \
    remove_host_reads/cow_host_bwa_hitsID.txt
#perl /hpf/largeprojects/ccmbio/nreinhardt/mugqic_pipelines/pipelines/metatranscriptomics/scripts/main_read_samout.pl cow host bwa pairs
#perl /hpf/largeprojects/ccmbio/nreinhardt/mugqic_pipelines/pipelines/metatranscriptomics/scripts/main_select_reads_fromfile.pl cow host bwa pairs

# cow1
/hpf/largeprojects/ccmbio/nreinhardt/mugqic_pipelines/pipelines/metatranscriptomics/scripts/split_reads_by_id.py \
    --fastq remove_rrna/cow1_qual_unique_n_rRNA.fastq \
    --id-file remove_host_reads/cow_host_bwa_IDs.txt \
    --included remove_host_reads/cow1_qual_unique_n_rRNA_host.fastq \
    --excluded remove_host_reads/cow1_qual_unique_n_rRNA_n_host.fastq
# cow2
/hpf/largeprojects/ccmbio/nreinhardt/mugqic_pipelines/pipelines/metatranscriptomics/scripts/split_reads_by_id.py \
    --fastq remove_rrna/cow2_qual_unique_n_rRNA.fastq \
    --id-file remove_host_reads/cow_host_bwa_IDs.txt \
    --included remove_host_reads/cow2_qual_unique_n_rRNA_host.fastq \
    --excluded remove_host_reads/cow2_qual_unique_n_rRNA_n_host.fastq
