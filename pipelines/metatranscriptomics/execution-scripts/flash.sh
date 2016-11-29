#!/bin/bash

mkdir flash

# TODO: make optional step
/hpf/largeprojects/ccmbio/nreinhardt/metatranscriptomics-test/reference-files/FLASH-1.2.11/flash -M 75 -p 64 -t 2 -d flash -o cow_qual trim/cow1_qual_paired.fastq trim/cow2_qual_paired.fastq
cat flash/cow_qual.extendedFrags.fastq flash/cow_qual.notCombined_1.fastq > flash/cow1_qual_all.fastq
cp flash/cow_qual.notCombined_2.fastq flash/cow2_qual_all.fastq
