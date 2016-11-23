#!/usr/bin/env bash

module load mugqic/bwa/0.7.12
module load mugqic/samtools/1.3
module load mugqic-pipelines/2.2.0
module load perl/5.20.1

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

#perl main_read_samout.pl cow assembly bwa pairs
perl /hpf/largeprojects/ccmbio/nreinhardt/mugqic_pipelines/pipelines/metatranscriptomics/scripts/main_read_samout.pl cow assembly bwa pairs

# cow1
python /hpf/largeprojects/ccmbio/nreinhardt/mugqic_pipelines/pipelines/metatranscriptomics/scripts/split_fastq_by_id.py \
    $MRNA_DIR/cow1_mRNA.fastq \
    $OUTPUT_DIR/cow_trinity_bwa_IDs.txt \
    $OUTPUT_DIR/cow1_mRNA_mappedreads.fastq \
    $OUTPUT_DIR/cow1_singletons.fastq
# cow2
python /hpf/largeprojects/ccmbio/nreinhardt/mugqic_pipelines/pipelines/metatranscriptomics/scripts/split_fastq_by_id.py \
    $MRNA_DIR/cow2_mRNA.fastq \
    $OUTPUT_DIR/cow_trinity_bwa_IDs.txt \
    $OUTPUT_DIR/cow2_mRNA_mappedreads.fastq \
    $OUTPUT_DIR/cow2_singletons.fastq
#perl main_select_reads_fromfile.pl cow assembly bwa pairs

# cow1
python /hpf/largeprojects/ccmbio/nreinhardt/mugqic_pipelines/pipelines/metatranscriptomics/scripts/main_get_sequence_length.py \
    --fastq $OUTPUT_DIR/cow1_singletons.fastq \
    --output $OUTPUT_DIR/cow1_singletons_IDs_length.txt
# cow2
python /hpf/largeprojects/ccmbio/nreinhardt/mugqic_pipelines/pipelines/metatranscriptomics/scripts/main_get_sequence_length.py \
    --fastq $OUTPUT_DIR/cow2_singletons.fastq \
    --output $OUTPUT_DIR/cow2_singletons_IDs_length.txt
#perl main_get_sequence_length.pl cow singletons

perl /hpf/largeprojects/ccmbio/nreinhardt/mugqic_pipelines/pipelines/metatranscriptomics/scripts/main_get_maptable_contig.pl \
    cow \
    assembly