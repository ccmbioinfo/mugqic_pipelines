#!/usr/bin/env bash

#PBS -N run-all
#PBS -l nodes=1:ppn=8
#PBS -l gres=localhd:20
#PBS -l vmem=50g
#PBS -l walltime=72:00:00

# Run all the steps

cd $PBS_O_WORKDIR

./format_fastq_headers_small.sh
./trim.sh
./flash.sh
./cluster_duplicates.sh
./remove_duplicates.sh
./remove_rrna.sh
./remove_host_reads.sh
./add_duplicates.sh
./trinity_assembly.sh
./index_contigs.sh
./map_reads.sh

./bwa_search_with_contigs.sh
./blat_search_with_contigs.sh
./diamond_search_with_contigs.sh

./bwa_search_with_singletons.sh
./blat_search_with_singletons.sh
#./diamond_search_with_singletons.sh
