#!/bin/bash
#PBS -l vmem=40g
#PBS -l nodes=1:ppn=6
#PBS -l walltime=48:00:00
#PBS -d /hpf/largeprojects/ccmbio/nreinhardt/metatranscriptomics-test/small-test

if [ "$PBS_ENVIRONMENT" = "PBS_BATCH" ]; then
    cd $PBS_O_WORKDIR
fi

module load infernal
module load perl

RFAM_PATH=/hpf/largeprojects/ccmbio/nreinhardt/metatranscriptomics-test/reference-files/Rfam.cm

mkdir remove_rrna
cmscan -o remove_rrna/cow1_rRNA.log --tblout remove_rrna/cow1_rRNA.infernalout --noali --notextw --rfam -E 0.001 $RFAM_PATH remove_duplicates/cow1_qual_all_unique.fasta
cmscan -o remove_rrna/cow2_rRNA.log --tblout remove_rrna/cow2_rRNA.infernalout --noali --notextw --rfam -E 0.001 $RFAM_PATH remove_duplicates/cow2_qual_all_unique.fasta


/hpf/largeprojects/ccmbio/nreinhardt/mugqic_pipelines/pipelines/metatranscriptomics/scripts/main_get_sequence_length.py
#perl main_get_sequence_length.pl cow rRNA
/hpf/largeprojects/ccmbio/nreinhardt/mugqic_pipelines/pipelines/metatranscriptomics/scripts/main_get_infernal_fromfile_1tophit.pl cow pairs 1 0.001 90
#perl main_get_infernal_fromfile_1tophit.pl cow pairs 1 0.001 90
/hpf/largeprojects/ccmbio/nreinhardt/mugqic_pipelines/pipelines/metatranscriptomics/scripts/main_select_reads_fromfile.pl
#perl main_select_reads_fromfile.pl cow rRNA infernal pairs
