#!/usr/bin/env bash

MRNA_DIR=add_duplicates
CONTIGS_DIR=index_contigs
OUTPUT_DIR=map_reads
mkdir $OUTPUT_DIR

bwa aln -t 4 $CONTIGS_DIR/cow_contigs.fasta $MRNA_DIR/cow1_mRNA.fastq > $OUTPUT_DIR/cow1_trinity.sai
bwa aln -t 4 $CONTIGS_DIR/cow_contigs.fasta $MRNA_DIR/cow2_mRNA.fastq > $OUTPUT_DIR/cow2_trinity.sai
bwa sampe $CONTIGS_DIR/cow_contigs.fasta $OUTPUT_DIR/cow1_trinity.sai $OUTPUT_DIR/cow2_trinity.sai  $MRNA_DIR/cow1_mRNA.fastq $MRNA_DIR/cow2_mRNA.fastq > $OUTPUT_DIR/cow_trinity.sam
samtools view -bS $OUTPUT_DIR/cow_trinity.sam | samtools sort -n -o $OUTPUT_DIR/cow_trinity.bam
samtools view $OUTPUT_DIR/cow_trinity.bam > $OUTPUT_DIR/cow_trinity.bwaout
#samtools view -F 4 $OUTPUT_DIR/cow_trinity.bam > $OUTPUT_DIR/cow_trinity.bwaout

perl /hpf/largeprojects/ccmbio/nreinhardt/mugqic_pipelines/pipelines/metatranscriptomics/scripts/main_read_samout.pl cow assembly bwa pairs
#perl main_read_samout.pl cow assembly bwa pairs
perl main_select_reads_fromfile.pl cow assembly bwa pairs
perl main_get_sequence_length.pl cow singletons

perl main_get_maptable_contig.pl cow assembly