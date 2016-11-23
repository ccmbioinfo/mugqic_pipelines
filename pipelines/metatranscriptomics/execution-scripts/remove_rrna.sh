#!/bin/bash
#PBS -l vmem=40g
#PBS -l nodes=1:ppn=6
#PBS -l walltime=48:00:00
#PBS -d /hpf/largeprojects/ccmbio/nreinhardt/mugqic_pipelines/pipelines/metatranscriptomics/execution-scripts

module load mugqic-pipelines/2.2.0
module load infernal
module load perl

RFAM_PATH=/hpf/largeprojects/ccmbio/nreinhardt/metatranscriptomics-test/reference-files/Rfam.cm

INPUT_DIR=remove_duplicates
OUTPUT_DIR=remove_rrna

mkdir $OUTPUT_DIR

cmscan -o $OUTPUT_DIR/cow1_rRNA.log --tblout $OUTPUT_DIR/cow1_rRNA.infernalout --noali --notextw --rfam -E 0.001 $RFAM_PATH %INPUT_DIR/cow1_qual_all_unique.fasta
cmscan -o $OUTPUT_DIR/cow2_rRNA.log --tblout $OUTPUT_DIR/cow2_rRNA.infernalout --noali --notextw --rfam -E 0.001 $RFAM_PATH %INPUT_DIR/cow2_qual_all_unique.fasta


# cow1
python /hpf/largeprojects/ccmbio/nreinhardt/mugqic_pipelines/pipelines/metatranscriptomics/scripts/main_get_sequence_length.py \
    --fasta $INPUT_DIR/cow1_qual_all_unique.fasta \
    --id-file $INPUT_DIR/cow1_qual_all_unique_IDs.txt \
    --output $OUTPUT_DIR/cow1_IDs_length.txt
# cow2
python /hpf/largeprojects/ccmbio/nreinhardt/mugqic_pipelines/pipelines/metatranscriptomics/scripts/main_get_sequence_length.py \
    --fasta $INPUT_DIR/cow2_qual_all_unique.fasta \
    --id-file $INPUT_DIR/cow2_qual_all_unique_IDs.txt \
    --output $OUTPUT_DIR/cow2_IDs_length.txt
#python /hpf/largeprojects/ccmbio/nreinhardt/mugqic_pipelines/pipelines/metatranscriptomics/scripts/main_get_sequence_length.py
#perl main_get_sequence_length.pl cow rRNA

perl /hpf/largeprojects/ccmbio/nreinhardt/mugqic_pipelines/pipelines/metatranscriptomics/scripts/main_get_infernal_fromfile_1tophit.pl cow pairs 1 0.001 90
#perl main_get_infernal_fromfile_1tophit.pl cow pairs 1 0.001 90

# cow1
python /hpf/largeprojects/ccmbio/nreinhardt/mugqic_pipelines/pipelines/metatranscriptomics/scripts/split_fastq_by_id.py \
    $INPUT_DIR/cow1_qual_all_unique.fastq \
    $OUTPUT_DIR/cow1_rRNA_infernal_IDs.txt \
    $OUTPUT_DIR/cow1_qual_unique_rRNA.fastq \
    $OUTPUT_DIR/cow1_qual_unique_n_rRNA.fastq
# cow2
python /hpf/largeprojects/ccmbio/nreinhardt/mugqic_pipelines/pipelines/metatranscriptomics/scripts/split_fastq_by_id.py \
    $INPUT_DIR/cow2_qual_all_unique.fastq \
    $OUTPUT_DIR/cow2_rRNA_infernal_IDs.txt \
    $OUTPUT_DIR/cow2_qual_unique_rRNA.fastq \
    $OUTPUT_DIR/cow2_qual_unique_n_rRNA.fastq
#perl /hpf/largeprojects/ccmbio/nreinhardt/mugqic_pipelines/pipelines/metatranscriptomics/scripts/main_select_reads_fromfile.pl
#perl main_select_reads_fromfile.pl cow rRNA infernal pairs
