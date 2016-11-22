#!/bin/bash

mkdir format_fastq_headers
/hpf/largeprojects/ccmbio/nreinhardt/mugqic_pipelines/pipelines/metatranscriptomics/scripts/main_add_subID_reads_fastq.pl \
    cow1_small.fastq \
    format_fastq_headers/cow1_new.fastq \
    cow2_small.fastq \
    format_fastq_headers/cow2_new.fastq
