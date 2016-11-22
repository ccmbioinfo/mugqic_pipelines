#!/usr/bin/env bash

#PBS -N run-all
#PBS -l nodes=1:ppn=8
#PBS -l gres=localhd:20
#PBS -l vmem=90g
#PBS -l walltime=100:00:00

cd $PBS_O_WORKDIR

./format_fastq_headers.sh
./flash.sh
./trim.sh
./cluster_duplicates.sh
./remoe_duplicates.sh
./remoe_rrna.sh
./remove_host_reads.sh
./add_duplicates.sh
./trinity_assembly.sh
./index_contigs.sh
./map_reads.sh