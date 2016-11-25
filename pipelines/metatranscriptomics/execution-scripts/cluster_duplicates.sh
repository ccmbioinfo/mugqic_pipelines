#!/bin/bash
set -x
module load seqtk
module load usearch
module load perl/5.20.1

mkdir cluster_duplicates


# Convert fastq to fasta
/hpf/largeprojects/ccmbio/nreinhardt/mugqic_pipelines/pipelines/metatranscriptomics/scripts/fastq_to_fasta.sh \
    flash/cow1_qual_all.fastq \
    cluster_duplicates/cow1_qual_all.fasta
/hpf/largeprojects/ccmbio/nreinhardt/mugqic_pipelines/pipelines/metatranscriptomics/scripts/fastq_to_fasta.sh \
    flash/cow2_qual_all.fastq \
    cluster_duplicates/cow2_qual_all.fasta
#seqtk seq -a flash/cow1_qual_all.fastq > cluster_duplicates/cow1_qual_all.fasta
#seqtk seq -a flash/cow2_qual_all.fastq > cluster_duplicates/cow2_qual_all.fasta


usearch --derep_fullseq --cluster cluster_duplicates/cow1_qual_all.fasta --seedsout cluster_duplicates/cow1_qual_all_unique.fasta --sizeout -uc cluster_duplicates/cow1_qual_all_unique.uc
usearch --derep_fullseq --cluster cluster_duplicates/cow2_qual_all.fasta --seedsout cluster_duplicates/cow2_qual_all_unique.fasta --sizeout -uc cluster_duplicates/cow2_qual_all_unique.uc

