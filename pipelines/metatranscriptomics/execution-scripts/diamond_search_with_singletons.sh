#!/usr/bin/env bash
# TODO

module load perl
module load blastp/2.2.27

blastdb=../reference-files
input_dir=blat_search_with_singletons
output_dir=diamond_search_with_singletons
mkdir $output_dir


# cow1
diamond blastx -p 8 -d $blastdb/nr.dmnd -q $input_dir/cow1_singletons_n_micro_cds_rest.fasta -a $output_dir/cow1_singletons_nr.matches -t $output_dir -e 10 -k 10
diamond view -a $output_dir/cow1_singletons_nr.matches.daa -o $output_dir/cow1_singletons_nr.diamondout -f tab
# cow2
diamond blastx -p 8 -d $blastdb/nr.dmnd -q $input_dir/cow2_singletons_n_micro_cds_rest.fasta -a $output_dir/cow2_singletons_nr.matches -t $output_dir -e 10 -k 10
diamond view -a $output_dir/cow2_singletons_nr.matches.daa -o $output_dir/cow2_singletons_nr.diamondout -f tab

# TODO no contigs
perl main_get_blast_fromfile_tophits.pl cow nr diamond contigs 1 100 85 65 60

perl main_sort_blastout_fromfile.pl cow nr diamond singletons 10

perl main_get_blast_fromfile_tophits.pl cow nr diamond singletons 1 100 85 65 60

# main_get_blast_fromfile_1topbachit.pl
perl /hpf/largeprojects/ccmbio/nreinhardt/mugqic_pipelines/pipelines/metatranscriptomics/scripts/main_get_blast_fromfile_1topbachit.pl \
    $output_dir/cow_singletons_nr_diamond_IDs.txt \
    $output_dir/cow_singletons_nr_diamond_pairs.txt \
    $output_dir/cow_singletons_nr_diamond_hitsID_sub.txt \
    $output_dir/cow_singletons_nr_diamond_pairs_sub.txt \
    $output_dir/cow_singletons_nr_diamond_hitsID_bacsub.txt
