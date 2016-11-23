#!/usr/bin/env bash

module load blat/35

blastdb=../reference-files
query_dir=bwa_search_with_contigs
output_dir=blat_search_with_contigs

blat -noHead -minIdentity=90 -minScore=50 $blastdb/microbial_all_cds.fasta $query_dir/cow_contigs_n_micro_cds.fasta -fine -q=rna -t=dna -out=blast8 $output_dir/cow_contigs.blatout

# TODO
perl main_sort_blastout_fromfile.pl cow n_micro_cds blat contigs 10
perl main_get_blast_fromfile_1tophit.pl cow micro_cds blat contigs 1 100 85 65 60
perl main_select_reads_fromfile.pl cow microgenes_blat blat contigs micro_cds
