#!/usr/bin/env bash

#PBS -l nodes=1:ppn=1
#PBS -l vmem=120g
#PBS -l walltime=23:00:00
#PBS -d /hpf/largeprojects/ccmbio/nreinhardt/mugqic_pipelines/pipelines/metatranscriptomics/execution-scripts
#PBS -joe /hpf/largeprojects/ccmbio/nreinhardt/mugqic_pipelines/pipelines/metatranscriptomics/execution-scripts

module load blat/35
module load mugqic-pipelines/2.2.0

blastdb=../reference-files
map_reads=map_reads
input_dir=bwa_search_with_contigs
output_dir=blat_search_with_contigs

mkdir $output_dir

# single fasta
#blat -noHead -minIdentity=90 -minScore=50 $blastdb/microbial_all_cds.fasta $input_dir/cow_contigs_n_micro_cds.fasta -fine -q=rna -t=dna -out=blast8 $output_dir/cow_contigs.blatout

# blat against 2 fastas
blat -noHead -minIdentity=90 -minScore=50 $blastdb/microbial_all_cds_1.fasta  $input_dir/cow_contigs_n_micro_cds.fasta -fine -q=rna -t=dna -out=blast8 $output_dir/cow_contigs_1.blatout
blat -noHead -minIdentity=90 -minScore=50 $blastdb/microbial_all_cds_2.fasta  $input_dir/cow_contigs_n_micro_cds.fasta -fine -q=rna -t=dna -out=blast8 $output_dir/cow_contigs_2.blatout
cat $output_dir/cow_contigs_1.blatout $output_dir/cow_contigs_2.blatout > $output_dir/cow_contigs_n_micro_cds.blatout

# main_sort_blastout_fromfile.pl
perl /hpf/largeprojects/ccmbio/nreinhardt/mugqic_pipelines/pipelines/metatranscriptomics/scripts/main_sort_blastout_fromfile.pl \
    $output_dir/cow_contigs_n_micro_cds.blatout \
    $output_dir/cow_contigs_n_micro_cds_sorted.blatout \
    10
#perl main_sort_blastout_fromfile.pl cow n_micro_cds blat contigs 10

# main_get_blast_fromfile_1tophit.pl
perl /hpf/largeprojects/ccmbio/nreinhardt/mugqic_pipelines/pipelines/metatranscriptomics/scripts/main_get_blast_fromfile_1tophit.pl \
    $output_dir/cow_contigs_n_micro_cds_sorted.blatout \
    $map_reads/cow_contigs_IDs_length.txt \
    $output_dir/cow_contigs_n_micro_cds_blat_IDs.txt \
    $output_dir/cow_contigs_n_micro_cds_blat_pairs.txt \
    $output_dir/cow_contigs_n_micro_cds_blat_hitsID.txt \
    1 100 85 65 60
#perl main_get_blast_fromfile_1tophit.pl cow micro_cds blat contigs 1 100 85 65 60

# split_reads_by_id.py
python /hpf/largeprojects/ccmbio/nreinhardt/mugqic_pipelines/pipelines/metatranscriptomics/scripts/split_reads_by_id.py \
    --fasta $input_dir/cow_contigs_n_micro_cds.fasta \
    --id-file $output_dir/cow_contigs_n_micro_cds_blat_IDs.txt \
    --included $output_dir/cow_contigs_n_micro_cds_blat.fasta \
    --excluded $output_dir/cow_contigs_n_micro_cds_rest.fasta
#perl main_select_reads_fromfile.pl cow microgenes_blat blat contigs micro_cds

