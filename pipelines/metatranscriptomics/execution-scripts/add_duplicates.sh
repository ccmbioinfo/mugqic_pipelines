#!/bin/bash

module load perl

id_dir=remove_duplicates
input_dir=remove_host_reads
output_dir=add_duplicates
mkdir $output_dir


# cow1
perl /hpf/largeprojects/ccmbio/nreinhardt/mugqic_pipelines/pipelines/metatranscriptomics/scripts/main_remove_derepli.pl \
    $id_dir/cow1_qual_all_unique_IDs.txt \
    $input_dir/cow1_qual_unique_n_rRNA_n_host.fastq \
    $output_dir/cow1_mRNA.fastq
# cow2
perl /hpf/largeprojects/ccmbio/nreinhardt/mugqic_pipelines/pipelines/metatranscriptomics/scripts/main_remove_derepli.pl \
    $id_dir/cow2_qual_all_unique_IDs.txt \
    $input_dir/cow2_qual_unique_n_rRNA_n_host.fastq \
    $output_dir/cow2_mRNA.fastq
