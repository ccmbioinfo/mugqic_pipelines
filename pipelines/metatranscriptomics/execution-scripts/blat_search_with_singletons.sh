#!/usr/bin/env bash

module load blat/35

blastdb=../reference-files
query_dir=bwa_search_with_singletons
output_dir=blat_search_with_singletons

# cow1
blat -noHead -minIdentity=90 -minScore=50 $blastdb/microbial_all_cds.fasta $query_dir/cow1_singletons_n_micro_cds.fasta -fine -q=rna -t=dna -out=blast8 cow1_singletons_n_micro_cds.blatout
# cow2
blat -noHead -minIdentity=90 -minScore=50 $blastdb/microbial_all_cds.fasta $query_dir/cow2_singletons_n_micro_cds.fasta -fine -q=rna -t=dna -out=blast8 cow2_singletons_n_micro_cds.blatout

# TODO
perl main_sort_blastout_fromfile.pl cow n_micro_cds blat singletons 10
perl main_get_blast_fromfile_1tophit.pl cow micro_cds blat singletons 1 100 85 65 60
perl main_select_reads_fromfile.pl cow microgenes_blat blat singletons micro_cds
