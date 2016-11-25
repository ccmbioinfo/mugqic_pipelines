#!/usr/bin/env bash

module load blat/35
module load mugqic-pipelines/2.2.0

blastdb=../reference-files
map_reads=map_reads
input_dir=bwa_search_with_singletons
output_dir=blat_search_with_singletons

mkdir $output_dir

# cow1
blat -noHead -minIdentity=90 -minScore=50 $blastdb/microbial_all_cds_1.fasta  $input_dir/cow1_singletons_n_micro_cds.fasta -fine -q=rna -t=dna -out=blast8 $output_dir/cow_singletons1_1.blatout
blat -noHead -minIdentity=90 -minScore=50 $blastdb/microbial_all_cds_2.fasta  $input_dir/cow1_singletons_n_micro_cds.fasta -fine -q=rna -t=dna -out=blast8 $output_dir/cow_singletons1_2.blatout
cat $output_dir/cow_singletons1_1.blatout $output_dir/cow_singletons1_2.blatout > $output_dir/cow1_singletons_n_micro_cds.blatout

# cow2
blat -noHead -minIdentity=90 -minScore=50 $blastdb/microbial_all_cds_1.fasta  $input_dir/cow2_singletons_n_micro_cds.fasta -fine -q=rna -t=dna -out=blast8 $output_dir/cow_singletons2_1.blatout
blat -noHead -minIdentity=90 -minScore=50 $blastdb/microbial_all_cds_2.fasta  $input_dir/cow2_singletons_n_micro_cds.fasta -fine -q=rna -t=dna -out=blast8 $output_dir/cow_singletons2_2.blatout
cat $output_dir/cow_singletons2_1.blatout $output_dir/cow_singletons2_2.blatout > $output_dir/cow2_singletons_n_micro_cds.blatout

cat $output_dir/cow1_singletons_n_micro_cds.blatout \
    $output_dir/cow2_singletons_n_micro_cds.blatout \
    > $output_dir/cow_singletons_n_micro_cds.blatout
perl /hpf/largeprojects/ccmbio/nreinhardt/mugqic_pipelines/pipelines/metatranscriptomics/scripts/main_sort_blastout_fromfile.pl \
    $output_dir/cow_singletons_n_micro_cds.blatout \
    $output_dir/cow_singletons_n_micro_cds_sorted.blatout \
    10
#perl main_sort_blastout_fromfile.pl cow n_micro_cds blat singletons 10

cat $map_reads/cow1_singletons_IDs_length.txt \
    $map_reads/cow1_singletons_IDs_length.txt \
    > $output_dir/cow_singletons_IDs_length.txt
perl /hpf/largeprojects/ccmbio/nreinhardt/mugqic_pipelines/pipelines/metatranscriptomics/scripts/main_get_blast_fromfile_1tophit.pl \
    $output_dir/cow_singletons_n_micro_cds_sorted.blatout \
    $output_dir/cow_singetons_IDs_length.txt \
    $output_dir/cow_singletons_n_micro_cds_blat_IDs.txt \
    $output_dir/cow_singletons_n_micro_cds_blat_pairs.txt \
    $output_dir/cow_singletons_n_micro_cds_blat_hitsID.txt \
    cow micro_cds blat singletons \
    1 100 85 65 60
#perl main_get_blast_fromfile_1tophit.pl cow micro_cds blat singletons 1 100 85 65 60

# cow1
python /hpf/largeprojects/ccmbio/nreinhardt/mugqic_pipelines/pipelines/metatranscriptomics/scripts/split_reads_by_id.py \
    --fasta $input_dir/cow2_singletons_n_micro_cds.fasta \
    --id-file $output_dir/cow_singletons_n_micro_cds_blat_IDs.txt \
    --included $output_dir/cow2_singletons_n_micro_cds_blat.fasta \
    --excluded $output_dir/cow2_singletons_n_micro_cds_rest.fasta
# cow2
python /hpf/largeprojects/ccmbio/nreinhardt/mugqic_pipelines/pipelines/metatranscriptomics/scripts/split_reads_by_id.py \
    --fasta $input_dir/cow2_singletons_n_micro_cds.fasta \
    --id-file $output_dir/cow_singletons_n_micro_cds_blat_IDs.txt \
    --included $output_dir/cow2_singletons_n_micro_cds_blat.fasta \
    --excluded $output_dir/cow2_singletons_n_micro_cds_rest.fasta
#perl main_select_reads_fromfile.pl cow microgenes_blat blat singletons micro_cds

