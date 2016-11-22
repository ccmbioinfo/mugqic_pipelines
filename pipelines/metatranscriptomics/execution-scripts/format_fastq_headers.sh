#!/bin/bash

module load perl

mkdir format_fastq_headers
/hpf/largeprojects/ccmbio/nreinhardt/mugqic_pipelines/pipelines/metatranscriptomics/scripts/main_add_subID_reads_fastq.pl \
    ../reference-files/cow1.fastq \
    format_fastq_headers/cow1_new.fastq \
    ../reference-files/cow2.fastq \
    format_fastq_headers/cow2_new.fastq
