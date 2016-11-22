#!/bin/bash

mkdir trim
java -jar /hpf/tools/centos6/trimmomatic/source/Trimmomatic-0.32/trimmomatic-0.32.jar \
    PE \
    format_fastq_headers/cow1_new.fastq format_fastq_headers/cow2_new.fastq trim/cow1_qual_paired.fastq trim/cow1_qual_unpaired.fastq trim/cow2_qual_paired.fastq trim/cow2_qual_unpaired.fastq \
    ILLUMINACLIP:/hpf/largeprojects/ccmbio/nreinhardt/metatranscriptomics-test/reference-files/UniVec_Core:2:30:10 \
    LEADING:3 \
    TRAILING:3 \
    SLIDINGWINDOW:4:15 \
    MINLEN:50
