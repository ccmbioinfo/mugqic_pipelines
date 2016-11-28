#!/usr/bin/env bash

module load mugqic-pipelines/2.2.0
module load mugqic/bwa/0.7.12
module load mugqic/samtools/1.3

INPUT_DIR=trinity_assembly
OUTPUT_DIR=index_contigs
mkdir $OUTPUT_DIR

cp $INPUT_DIR/cow_contigs.fasta $OUTPUT_DIR/cow_contigs.fasta

bwa index -a bwtsw $OUTPUT_DIR/cow_contigs.fasta
samtools faidx $OUTPUT_DIR/cow_contigs.fasta
