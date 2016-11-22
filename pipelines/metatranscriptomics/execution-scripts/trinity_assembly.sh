#!/bin/bash -x

#PBS -N Job_trinity
#PBS -l nodes=1:ppn=8
#PBS -l gres=localhd:20
#PBS -l vmem=90g
#PBS -l walltime=48:00:00
#PBS -joe /hpf/largeprojects/ccmbio/nreinhardt/mugqic_pipelines/pipelines/metatranscriptomics/execution-scripts


module load trinityrnaseq/2.1.1
module load bowtie/1.0.1
module load samtools/1.2

my_WORKDIR=/hpf/largeprojects/ccmbio/nreinhardt/mugqic_pipelines/pipelines/metatranscriptomics/execution-scripts
cd ${my_WORKDIR}

indir=add_duplicates
outdir=trinity_assembly

rm -fr trinity_assembly
mkdir trinity_assembly

Trinity --seqType fq --left $indir/cow1_mRNA.fastq --right $indir/cow2_mRNA.fastq --CPU 8 --max_memory 10G --min_contig_length 75 --output ${my_WORKDIR}/${outdir} --no_version_check
#Trinity --seqType fq --left $indir/cow1_mRNA.fastq --right $indir/cow2_mRNA.fastq --CPU 8 --max_memory 10G --min_contig_length 75 --output ${my_WORKDIR}/${outdir} --full_cleanup --no_version_check

cp ${my_WORKDIR}/${outdir}/Trinity.fasta ${my_WORKDIR}/$outdir/cow_contigs.fasta
