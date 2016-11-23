#!/usr/bin/env bash

module load mugqic-pipelines/2.2.0
module load perl
module load mugqic/bwa/0.7.12
module load mugqic/samtools/1.3
module load blastp/2.2.27

input_dir=map_reads
contigs_dir=index_contigs
output_dir=bwa_search_with_contigs
reference_dir=../reference-files

mkdir $output_dir


# Contig searches
bwa aln -t 4 $reference_dir/microbial_all_cds.fasta $contigs_dir/cow_contigs.fasta > $output_dir/cow_contigs.sai
bwa samse $reference_dir/microbial_all_cds.fasta $output_dir/cow_contigs.sai $contigs_dir/cow_contigs.fasta > $output_dir/cow_contigs.sam
samtools view -bS $output_dir/cow_contigs.sam | samtools sort -n -o $output_dir/cow_contigs.bam
samtools view -F 4 $output_dir/cow_contigs.bam > $output_dir/cow_contigs_micro_cds.bwaout

perl /hpf/largeprojects/ccmbio/nreinhardt/mugqic_pipelines/pipelines/metatranscriptomics/scripts/main_read_samout.pl \
    $output_dir/cow_contigs_micro_cds.bwaout \
    $output_dir/cow_contigs_micro_cds_bwa_IDs.txt \
    $output_dir/cow_contigs_micro_cds_bwa_pairs.txt \
    $output_dir/cow_contigs_micro_cds_bwa_hitsID.txt
#perl main_read_samout.pl cow microgenes bwa contigs micro_cds

python /hpf/largeprojects/ccmbio/nreinhardt/mugqic_pipelines/pipelines/metatranscriptomics/scripts/split_reads_by_id.py \
    --fasta $contigs_dir/cow_contigs.fasta \
    --id-file $output_dir/cow_contigs_micro_cds_bwa_IDs.txt \
    --included $output_dir/cow_contigs_micro_cds.fasta \
    --excluded $output_dir/cow_contigs_n_micro_cds.fasta
#perl main_select_reads_fromfile.pl cow microgenes bwa contigs micro_cds


