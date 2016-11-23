#!/usr/bin/env bash

module load mugqic/bwa/0.7.12
module load mugqic/samtools/1.3
module load mugqic-pipelines/2.2.0
module load perl/5.20.1

mrna_dir=add_duplicates
contigs_dir=index_contigs
output_dir=map_reads
mkdir $output_dir

bwa aln -t 4 $contigs_dir/cow_contigs.fasta $mrna_dir/cow1_mRNA.fastq > $output_dir/cow1_trinity.sai
bwa aln -t 4 $contigs_dir/cow_contigs.fasta $mrna_dir/cow2_mRNA.fastq > $output_dir/cow2_trinity.sai
bwa sampe $contigs_dir/cow_contigs.fasta $output_dir/cow1_trinity.sai $output_dir/cow2_trinity.sai  $mrna_dir/cow1_mRNA.fastq $mrna_dir/cow2_mRNA.fastq > $output_dir/cow_trinity.sam
samtools view -bS $output_dir/cow_trinity.sam | samtools sort -n -o $output_dir/cow_trinity.bam
samtools view $output_dir/cow_trinity.bam > $output_dir/cow_trinity.bwaout
#samtools view -F 4 $OUTPUT_DIR/cow_trinity.bam > $OUTPUT_DIR/cow_trinity.bwaout

#perl main_read_samout.pl cow assembly bwa pairs
#perl /hpf/largeprojects/ccmbio/nreinhardt/mugqic_pipelines/pipelines/metatranscriptomics/scripts/main_read_samout.pl cow assembly bwa pairs
perl /hpf/largeprojects/ccmbio/nreinhardt/mugqic_pipelines/pipelines/metatranscriptomics/scripts/main_read_samout.pl \
    $output_dir/cow_trinity.bwaout \
    $output_dir/cow_trinity_bwa_IDs.txt \
    $output_dir/cow_trinity_bwa_pairs.txt \
    $output_dir/cow_trinity_bwa_hitsID.txt


# cow1
python /hpf/largeprojects/ccmbio/nreinhardt/mugqic_pipelines/pipelines/metatranscriptomics/scripts/split_reads_by_id.py \
    --fastq $mrna_dir/cow1_mRNA.fastq \
    --id-file $output_dir/cow_trinity_bwa_IDs.txt \
    --included $output_dir/cow1_mRNA_mappedreads.fastq \
    --excluded $output_dir/cow1_singletons.fastq
# cow2
python /hpf/largeprojects/ccmbio/nreinhardt/mugqic_pipelines/pipelines/metatranscriptomics/scripts/split_reads_by_id.py \
    --fastq $mrna_dir/cow2_mRNA.fastq \
    --id-file $output_dir/cow_trinity_bwa_IDs.txt \
    --included $output_dir/cow2_mRNA_mappedreads.fastq \
    --excluded $output_dir/cow2_singletons.fastq
#perl main_select_reads_fromfile.pl cow assembly bwa pairs

# cow1
python /hpf/largeprojects/ccmbio/nreinhardt/mugqic_pipelines/pipelines/metatranscriptomics/scripts/main_get_sequence_length.py \
    --fastq $output_dir/cow1_singletons.fastq \
    --output $output_dir/cow1_singletons_IDs_length.txt
# cow2
python /hpf/largeprojects/ccmbio/nreinhardt/mugqic_pipelines/pipelines/metatranscriptomics/scripts/main_get_sequence_length.py \
    --fastq $output_dir/cow2_singletons.fastq \
    --output $output_dir/cow2_singletons_IDs_length.txt
#perl main_get_sequence_length.pl cow singletons

perl /hpf/largeprojects/ccmbio/nreinhardt/mugqic_pipelines/pipelines/metatranscriptomics/scripts/main_get_maptable_contig.pl \
    cow \
    assembly