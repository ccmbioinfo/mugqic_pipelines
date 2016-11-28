#!/usr/bin/env bash

#PBS -l nodes=1:ppn=8
#PBS -l vmem=60g
#PBS -l walltime=23:00:00
#PBS -d /hpf/largeprojects/ccmbio/nreinhardt/mugqic_pipelines/pipelines/metatranscriptomics/execution-scripts
#PBS -joe /hpf/largeprojects/ccmbio/nreinhardt/mugqic_pipelines/pipelines/metatranscriptomics/execution-scripts

module load perl
module load blastp/2.2.27

blastdb=../reference-files
map_reads=map_reads
input_dir=blat_search_with_contigs
output_dir=diamond_search_with_contigs
mkdir $output_dir

diamond blastx -p 8 -d $blastdb/nr -q $input_dir/cow_contigs_n_micro_cds_rest.fasta -a $output_dir/cow_contigs_nr.matches -t $output_dir -e 10 -k 10
diamond view -a $output_dir/cow_contigs_nr.matches.daa -o $output_dir/cow_contigs_nr.diamondout -f tab

# main_get_blast_fromfile_tophits.pl
# TODO: pass in read_type, etc.
perl /hpf/largeprojects/ccmbio/nreinhardt/mugqic_pipelines/pipelines/metatranscriptomics/scripts/main_get_blast_fromfile_tophits.pl \
    $output_dir/cow_contigs_nr.diamondout \
    $map_reads/cow_contigs_IDs_length.txt \
    $output_dir/cow_contigs_nr_diamond_IDs.txt \
    $output_dir/cow_contigs_nr_diamond_pairs.txt \
    $output_dir/cow_contigs_nr_diamond_hitsID.txt \
    1 100 85 65 60
#perl main_get_blast_fromfile_tophits.pl cow nr diamond contigs 1 100 85 65 60

# main_get_blast_fromfile_1topbachit.pl
perl /hpf/largeprojects/ccmbio/nreinhardt/mugqic_pipelines/pipelines/metatranscriptomics/scripts/main_get_blast_fromfile_1topbachit.pl \
    $output_dir/cow_contigs_nr_diamond_IDs.txt \
    $output_dir/cow_contigs_nr_diamond_pairs.txt \
    $output_dir/cow_contigs_nr_diamond_hitsID_sub.txt \
    $output_dir/cow_contigs_nr_diamond_pairs_sub.txt \
    $output_dir/cow_contigs_nr_diamond_hitsID_bacsub.txt
#perl main_get_blast_fromfile_1topbachit.pl cow nr diamond contigs

