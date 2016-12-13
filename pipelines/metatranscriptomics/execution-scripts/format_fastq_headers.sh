#!/bin/bash
# $1 - path to fastq1
# $2 - path to fastq2

module load perl

mkdir format_fastq_headers
/hpf/largeprojects/ccmbio/nreinhardt/mugqic_pipelines/pipelines/metatranscriptomics/scripts/main_add_subID_reads_fastq.pl \
    $1 \
    format_fastq_headers/cow1_new.fastq \
    $2 \
    format_fastq_headers/cow2_new.fastq
