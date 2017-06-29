#!/bin/bash
# Exit immediately on error
set -eu -o pipefail

#-------------------------------------------------------------------------------
# Metatranscriptomics PBSScheduler Job Submission Bash script
# Version: 2.2.1-beta
# Created on: 2017-06-28T17:17:54
# Steps:
#   format_fastq_headers: 0 job... skipping
#   trimmomatic: 0 job... skipping
#   merge_overlapping_reads: 0 job... skipping
#   fastq_to_fasta: 0 job... skipping
#   cluster_duplicates: 0 job... skipping
#   remove_duplicates: 0 job... skipping
#   cmscan: 0 job... skipping
#   identify_rrna: 0 job... skipping
#   remove_rrna: 0 job... skipping
#   align_to_host: 0 job... skipping
#   identify_host_reads: 0 job... skipping
#   remove_host_reads: 0 job... skipping
#   return_duplicates: 0 job... skipping
#   trinity: 0 job... skipping
#   index_contigs: 0 job... skipping
#   align_to_contigs: 0 job... skipping
#   identify_contigs_reads: 0 job... skipping
#   extract_singletons: 0 job... skipping
#   get_mapping_table: 0 job... skipping
#   bwa_align_contigs: 0 job... skipping
#   bwa_identify_contigs: 0 job... skipping
#   bwa_contigs_select_reads: 0 job... skipping
#   bwa_align_singletons: 0 job... skipping
#   bwa_identify_singletons: 0 job... skipping
#   bwa_singletons_select_reads: 0 job... skipping
#   blat_search_contigs: 0 job... skipping
#   blat_search_singletons: 0 job... skipping
#   process_contigs: 0 job... skipping
#   process_singletons: 0 job... skipping
#   diamond_align_contigs: 0 job... skipping
#   diamond_align_singletons: 0 job... skipping
#   diamond_contigs_get_tophits: 0 job... skipping
#   diamond_singletons_get_tophits: 0 job... skipping
#   generate_microbial_sequence: 0 job... skipping
#   get_topbachit_contigs: 0 job... skipping
#   get_topbachit_singletons: 0 job... skipping
#   generate_nr_sequence: 0 job... skipping
#   align_genes_ecoli: 0 job... skipping
#   align_proteins_ecoli: 0 job... skipping
#   combine_ppi_results: 0 job... skipping
#   get_taxID_microbial: 0 job... skipping
#   get_taxID_nr: 0 job... skipping
#   get_phylum_microbial: 0 job... skipping
#   get_phylum_nr: 0 job... skipping
#   get_mapped_geneIDs_microbial: 0 job... skipping
#   get_mapped_geneIDs_nr: 2 jobs
#   get_mapped_gene_table_microbial: 0 job... skipping
#   get_mapped_gene_table_nr: 1 job
#   generate_RPKM: 2 jobs
#   TOTAL: 5 jobs
#-------------------------------------------------------------------------------

OUTPUT_DIR=/hpf/projects/brudno/daniel/mugqic_pipelines/pipelines/metatranscriptomics/output
JOB_OUTPUT_DIR=$OUTPUT_DIR/job_output
TIMESTAMP=`date +%FT%H.%M.%S`
JOB_LIST=$JOB_OUTPUT_DIR/Metatranscriptomics_job_list_$TIMESTAMP
mkdir -p $OUTPUT_DIR
cd $OUTPUT_DIR


#-------------------------------------------------------------------------------
# STEP: get_mapped_geneIDs_nr
#-------------------------------------------------------------------------------
STEP=get_mapped_geneIDs_nr
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: get_mapped_geneIDs_nr_1_JOB_ID: get_mapped_geneIDs_nr.cow_readset.combine_pairs
#-------------------------------------------------------------------------------
JOB_NAME=get_mapped_geneIDs_nr.cow_readset.combine_pairs
JOB_DEPENDENCIES=
JOB_DONE=job_output/get_mapped_geneIDs_nr/get_mapped_geneIDs_nr.cow_readset.combine_pairs.05e213eeb947c67a0c747ea4e1bb9765.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'get_mapped_geneIDs_nr.cow_readset.combine_pairs.05e213eeb947c67a0c747ea4e1bb9765.mugqic.done'
TMPDIR=/localhd/$PBS_JOBID
module load mugqic-pipelines/2.2.0
cat contigs/cow_readset/cow_readset.contigs.nr_diamond_pairs_sub.txt contigs/cow_readset/cow_readset.singletons.nr_diamond_pairs_sub.txt > contigs/cow_readset/cow_readset.nr_all_sub_combined_pairs_sub.txt
get_mapped_geneIDs_nr.cow_readset.combine_pairs.05e213eeb947c67a0c747ea4e1bb9765.mugqic.done
)
get_mapped_geneIDs_nr_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -W umask=0002 -W group_list=ccm -l vmem=10g,mem=10g -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=12:00:00  -l nodes=1:ppn=1 | grep "[0-9]")
usleep 500
echo "$get_mapped_geneIDs_nr_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: get_mapped_geneIDs_nr_2_JOB_ID: get_mapped_geneIDs_nr.cow_readset.map_geneIDs
#-------------------------------------------------------------------------------
JOB_NAME=get_mapped_geneIDs_nr.cow_readset.map_geneIDs
JOB_DEPENDENCIES=$get_mapped_geneIDs_nr_1_JOB_ID
JOB_DONE=job_output/get_mapped_geneIDs_nr/get_mapped_geneIDs_nr.cow_readset.map_geneIDs.ff8fed1a243ac0efe9472d53f2d4dec8.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'get_mapped_geneIDs_nr.cow_readset.map_geneIDs.ff8fed1a243ac0efe9472d53f2d4dec8.mugqic.done'
TMPDIR=/localhd/$PBS_JOBID
module load mugqic-pipelines/2.2.0
perl /hpf/projects/brudno/daniel/mugqic_pipelines/pipelines/metatranscriptomics/scripts/main_get_mapped_genesID_counts.pl contigs/cow_readset/cow_readset.nr_all_sub_IDs_length.txt contigs/cow_readset/cow_readset.contigs.IDs_length.txt contigs/cow_readset/cow_readset.nr_all_sub_combined_pairs_sub.txt contigs/cow_readset/cow_readset.nr_all_sub_IDs_counts.txt
get_mapped_geneIDs_nr.cow_readset.map_geneIDs.ff8fed1a243ac0efe9472d53f2d4dec8.mugqic.done
)
get_mapped_geneIDs_nr_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -W umask=0002 -W group_list=ccm -l vmem=10g,mem=10g -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=12:00:00  -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
usleep 500
echo "$get_mapped_geneIDs_nr_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: get_mapped_gene_table_nr
#-------------------------------------------------------------------------------
STEP=get_mapped_gene_table_nr
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: get_mapped_gene_table_nr_1_JOB_ID: get_mapped_gene_table_nr.cow_readset.get_mapped_gene_table_nr
#-------------------------------------------------------------------------------
JOB_NAME=get_mapped_gene_table_nr.cow_readset.get_mapped_gene_table_nr
JOB_DEPENDENCIES=$get_mapped_geneIDs_nr_2_JOB_ID
JOB_DONE=job_output/get_mapped_gene_table_nr/get_mapped_gene_table_nr.cow_readset.get_mapped_gene_table_nr.6770928f8982b5ff8ad4edc0f10e8c44.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'get_mapped_gene_table_nr.cow_readset.get_mapped_gene_table_nr.6770928f8982b5ff8ad4edc0f10e8c44.mugqic.done'
TMPDIR=/localhd/$PBS_JOBID
module load mugqic-pipelines/2.2.0
perl /hpf/projects/brudno/daniel/mugqic_pipelines/pipelines/metatranscriptomics/scripts/main_get_mapped_gene_table.pl contigs/cow_readset/cow_readset.nr_all_sub_IDs_length.txt contigs/cow_readset/cow_readset.nr_all_sub_IDs_map_taxid_phylum.txt contigs/cow_readset/cow_readset.nr_all_sub_IDs_counts.txt contigs/cow_readset/cow_readset.PPI_pairs.txt contigs/cow_readset/cow_readset.nr_all_sub_table_counts.txt
get_mapped_gene_table_nr.cow_readset.get_mapped_gene_table_nr.6770928f8982b5ff8ad4edc0f10e8c44.mugqic.done
)
get_mapped_gene_table_nr_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -W umask=0002 -W group_list=ccm -l vmem=10g,mem=10g -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=12:00:00  -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
usleep 500
echo "$get_mapped_gene_table_nr_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# STEP: generate_RPKM
#-------------------------------------------------------------------------------
STEP=generate_RPKM
mkdir -p $JOB_OUTPUT_DIR/$STEP


#-------------------------------------------------------------------------------
# JOB: generate_RPKM_1_JOB_ID: generate_RPKM.cow_readset.combine_counts
#-------------------------------------------------------------------------------
JOB_NAME=generate_RPKM.cow_readset.combine_counts
JOB_DEPENDENCIES=$get_mapped_gene_table_nr_1_JOB_ID
JOB_DONE=job_output/generate_RPKM/generate_RPKM.cow_readset.combine_counts.4856c52be37604708fb5df797e042ad5.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'generate_RPKM.cow_readset.combine_counts.4856c52be37604708fb5df797e042ad5.mugqic.done'
TMPDIR=/localhd/$PBS_JOBID
module load mugqic-pipelines/2.2.0
cat contigs/cow_readset/cow_readset.microbial_cds_sub_table_counts.txt contigs/cow_readset/cow_readset.nr_all_sub_table_counts.txt > contigs/cow_readset/cow_readset.table_counts_all
generate_RPKM.cow_readset.combine_counts.4856c52be37604708fb5df797e042ad5.mugqic.done
)
generate_RPKM_1_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -W umask=0002 -W group_list=ccm -l vmem=10g,mem=10g -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=12:00:00  -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
usleep 500
echo "$generate_RPKM_1_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# JOB: generate_RPKM_2_JOB_ID: generate_RPKM.cow_readset.generate_RPKM
#-------------------------------------------------------------------------------
JOB_NAME=generate_RPKM.cow_readset.generate_RPKM
JOB_DEPENDENCIES=$generate_RPKM_1_JOB_ID
JOB_DONE=job_output/generate_RPKM/generate_RPKM.cow_readset.generate_RPKM.17b9c322cae451dc723c42b89449de4c.mugqic.done
JOB_OUTPUT_RELATIVE_PATH=$STEP/${JOB_NAME}_$TIMESTAMP.o
JOB_OUTPUT=$JOB_OUTPUT_DIR/$JOB_OUTPUT_RELATIVE_PATH
COMMAND=$(cat << 'generate_RPKM.cow_readset.generate_RPKM.17b9c322cae451dc723c42b89449de4c.mugqic.done'
TMPDIR=/localhd/$PBS_JOBID
module load mugqic-pipelines/2.2.0
perl /hpf/projects/brudno/daniel/mugqic_pipelines/pipelines/metatranscriptomics/scripts/main_get_mapped_gene_table_RPKM.pl contigs/cow_readset/cow_readset.table_counts_all contigs/cow_readset/cow_readset.table_RPKM_all.txt
generate_RPKM.cow_readset.generate_RPKM.17b9c322cae451dc723c42b89449de4c.mugqic.done
)
generate_RPKM_2_JOB_ID=$(echo "rm -f $JOB_DONE && $COMMAND
MUGQIC_STATE=\$PIPESTATUS
echo MUGQICexitStatus:\$MUGQIC_STATE
if [ \$MUGQIC_STATE -eq 0 ] ; then touch $JOB_DONE ; fi
exit \$MUGQIC_STATE" | \
qsub -m ae -W umask=0002 -W group_list=ccm -l vmem=10g,mem=10g -d $OUTPUT_DIR -j oe -o $JOB_OUTPUT -N $JOB_NAME -l walltime=12:00:00  -l nodes=1:ppn=1 -W depend=afterok:$JOB_DEPENDENCIES | grep "[0-9]")
usleep 500
echo "$generate_RPKM_2_JOB_ID	$JOB_NAME	$JOB_DEPENDENCIES	$JOB_OUTPUT_RELATIVE_PATH" >> $JOB_LIST


#-------------------------------------------------------------------------------
# Call home with pipeline statistics
#-------------------------------------------------------------------------------
wget "http://mugqic.hpc.mcgill.ca/cgi-bin/pipeline.cgi?hostname=qlogin7&ip=10.10.2.2&pipeline=Metatranscriptomics&steps=format_fastq_headers,trimmomatic,merge_overlapping_reads,fastq_to_fasta,cluster_duplicates,remove_duplicates,cmscan,identify_rrna,remove_rrna,align_to_host,identify_host_reads,remove_host_reads,return_duplicates,trinity,index_contigs,align_to_contigs,identify_contigs_reads,extract_singletons,get_mapping_table,bwa_align_contigs,bwa_identify_contigs,bwa_contigs_select_reads,bwa_align_singletons,bwa_identify_singletons,bwa_singletons_select_reads,blat_search_contigs,blat_search_singletons,process_contigs,process_singletons,diamond_align_contigs,diamond_align_singletons,diamond_contigs_get_tophits,diamond_singletons_get_tophits,generate_microbial_sequence,get_topbachit_contigs,get_topbachit_singletons,generate_nr_sequence,align_genes_ecoli,align_proteins_ecoli,combine_ppi_results,get_taxID_microbial,get_taxID_nr,get_phylum_microbial,get_phylum_nr,get_mapped_geneIDs_microbial,get_mapped_geneIDs_nr,get_mapped_gene_table_microbial,get_mapped_gene_table_nr,generate_RPKM&samples=1" --quiet --output-document=/dev/null

