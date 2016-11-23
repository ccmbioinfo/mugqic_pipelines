#!/usr/bin/env bash

#PBS -N run-all
#PBS -l nodes=1:ppn=2
#PBS -l gres=localhd:20
#PBS -l vmem=40g
#PBS -l walltime=23:00:00

# Run all the steps

cd $PBS_O_WORKDIR

./format_fastq_headers_small.sh
./flash.sh
./trim.sh
./cluster_duplicates.sh
./remove_duplicates.sh
./remove_rrna.sh
./remove_host_reads.sh
./add_duplicates.sh
./trinity_assembly.sh
./index_contigs.sh
#./map_reads.sh