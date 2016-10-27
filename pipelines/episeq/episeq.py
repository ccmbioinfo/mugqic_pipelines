#!/usr/bin/env python
"""
Epigenetics pipeline for RRBS/WGBS data
=======================================

Background
----------
EpiSeq is a differential analysis pipeline for BS-seq sequencing. Currently only RRBS and WGBS datasets are
tested to work with this pipeline. Similar to the other MUGQIC pipeline series, EpiSeq uses two metadata files to
set up the pipeline. The design file is used to group samples into case vs control. The readsets files refer to
which input files correspond to each same and, if applicable, which pairs of input files correspond to paired reads.
A readset is considered to be one instance/lane/segment of sequencing data. This is often used when multiple libraries
are used for a given sample, or if multiplexing was done. These techniques tend to generate multiple sets of data
for a given sample. The readset file allows users to specify these relationships.

This simple pipeline uses Bismark, Trim Galore!, and Picard to process data. See the README file for more information
about the steps in the pipeline.

Implementation
--------------
The pipeline is arranged to gradually decrease the number of input files when it make sense to do so. The trimming
and QC step would be too slow if we merge readsets together, so we decided to trim each file individually. Then,
the alignment phase allows us to convert out fastQ files to BAM files. This allows us to merge paired datasets into
one BAM file per readset without explicitly doing so. After this, readsets are combined to run the differential
analysis steps.

Steps (1 and 2) and (6 and 7) run concurrently while all other steps run sequentially.
"""

# Python Standard Modules
import argparse
import collections
import logging
import os
import re
import sys

# Append mugqic_pipelines directory to Python library path
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0])))))

# MUGQIC Modules
from core.config import *
from core.job import *
from core.pipeline import *
from bfx.readset import *
from bfx.design import *

from bfx import rmarkdown
from pipelines import common
import utils

log = logging.getLogger(__name__)


class Episeq(common.Illumina):
    """
    The Episeq pipeline takes FASTQ or BAM files (unsorted) as input as well as two metadata files and a configuration
    file. Refer to the user guide for more information on running the pipeline.
    """

    def __init__(self):
        self.argparser.add_argument("-d", "--design", help="design file", type=file)
        super(Episeq, self).__init__()

    def merge_fastq(self):

        """
        This step merges FASTQ files for samples with multiple readsets or creates symbolic links of FASTQ
        files for samples with only a single readset
        """

        jobs = []
        for sample in self.samples:
            merge_directory = os.path.join("merged", sample.name)
            merge_prefix = os.path.join(merge_directory, sample.name)
            mkdir_job = Job(command="mkdir -p " + merge_directory)
            # If one readset in a sample is of a particular type,
            # then all readsets are assumed to be of that type as well
            run_type = sample.readsets[0].run_type

            job = []
            # Samples with only one readset do not require to be merged.
            # A symbolic link to the FASTQ files are made instead
            if len(sample.readsets) == 1:
                target_fastq1 = sample.readsets[0].fastq1
                if run_type == "PAIRED_END":
                    target_fastq2 = sample.readsets[0].fastq2
                elif run_type == "SINGLE_END":
                    target_fastq2 = ""

                job = concat_jobs([
                    mkdir_job,
                    Job([target_fastq1], [merge_prefix + "_R1.fastq.gz"],
                        command="ln -s -f " + target_fastq1 + " " + merge_prefix + "_R1.fastq.gz"),
                    Job([target_fastq2], [merge_prefix + "_R2.fastq.gz" if run_type == "PAIRED_END" else ""],
                        command="ln -s -f " + target_fastq2 + " " + merge_prefix + "_R2.fastq.gz")
                ], name="symlink_fastq." + sample.name)

            # Samples with more than one readsets will be merged by their forward and/or reverse FASTQ files
            elif len(sample.readsets) > 1:
                fastq1 = [readset.fastq1 for readset in sample.readsets]
                if run_type == "PAIRED_END":
                    fastq2 = [readset.fastq2 for readset in sample.readsets]
                elif run_type == "SINGLE_END":
                    fastq2 = ""

                job = concat_jobs([
                    mkdir_job,
                    Job([fastq for fastq in fastq1],
                        [merge_prefix + "_R1.fastq.gz"],
                        command="cat {fastq1} > {merged_fastq1}".format(
                            fastq1=" ".join(fastq1),
                            merged_fastq1=merge_prefix + "_R1.fastq.gz")),
                    Job([fastq for fastq in fastq2],
                        [merge_prefix + "_R2.fastq.gz" if run_type == "PAIRED_END" else ""],
                        command="cat {fastq2} > {merged_fastq2}".format(
                            fastq2=" ".join(fastq2),
                            merged_fastq2=merge_prefix + "_R2.fastq.gz"))
                ], name="merge_fastq." + sample.name)

            jobs.append(job)

        return jobs

    def trim_galore(self):
        """
        This step trims raw FASTQ files for quality control using Trim Galore!
        This is a pre-proccessing step to ensure quality control.

        :return jobs: A list of jobs that needs to be executed in this step.
        :rtype list(Job):
        """

        jobs = []
        for sample in self.samples:
            merge_prefix = os.path.join("merged", sample.name, sample.name)
            trim_directory = os.path.join("trimmed", sample.name)
            trim_prefix = os.path.join(trim_directory, sample.name)
            run_type = sample.readsets[0].run_type
            protocol = sample.readsets[0].library
            output_files = []

            # Trim Galoree has no built in option to change the filenames of the output
            # Below are the default output names when running in paired or single mode
            if run_type == "PAIRED_END":
                input_files = [merge_prefix + "_R1.fastq.gz",
                               merge_prefix + "_R2.fastq.gz"]
                output_files = [trim_prefix + "_R1_val_1.fq.gz",
                                trim_prefix + "_R2_val_2.fq.gz"]
            elif run_type == "SINGLE_END":
                input_files = [merge_prefix + "_R1.fastq.gz"]
                output_files = [trim_prefix + "_R1_trimmed.fq.gz"]

            mkdir_job = Job(command="mkdir -p " + trim_directory)
            job = concat_jobs([
                mkdir_job,
                Job(
                    input_files,
                    output_files,
                    command="""\
module load trim_galore/0.4.1
trim_galore {protocol} {library_type} {non_directional} {other_options} --output_dir {directory} {fastq1} {fastq2}
    """.format(
                        library_type="--paired" if run_type == "PAIRED_END" else "",
                        protocol='--rrbs' if protocol == 'RRBS' else '',
                        non_directional='--non_directional' if protocol == 'RRBS' and run_type == 'PAIRED_END' else '',
                        other_options=config.param("trim_galore", "other_options"),
                        directory=trim_directory,
                        fastq1=input_files[0],
                        fastq2=input_files[1] if run_type == "PAIRED_END" else ""
                    )  # Add --fastqc?
                )], name="trim_galore." + sample.name)
            jobs.append(job)

        return jobs

    def bismark_prepare_genome(self):
        """
        Bismark requires a processed reference genome to compare with the epigenome. This step can take several hours,
        depending on the size of the genome. It creates an index of bisulfite conversions and takes make space than
        the genome itself.

        The genome should be in fasta format and note that this step always copies the genome to the output directory.

        Input: A reference sequence file as specified by user. Configuration is set in the episeq.ini file.
        Output: A directory called Bisulfite_Genome.

        :return jobs: A list of jobs that needs to be executed in this step.
        :rtype list(Job):
        """
        ref_seq = config.param('bismark_prepare_genome', 'genome_file', type='filepath')
        local_ref_seq = os.path.join('bismark_prepare_genome', os.path.basename(ref_seq))
        output_idx = "Bisulfite_Genome"

        main_job = Job([local_ref_seq], [output_idx],
                       [['bismark_prepare_genome', 'module_bowtie2'],
                        ['bismark_prepare_genome', 'module_samtools']],
                       command="module load bismark/0.15; bismark_genome_preparation --verbose --bowtie2 .",
                       name="bismark_prepare_genome")

        if not os.path.isfile(local_ref_seq):
            link_job = Job([ref_seq], [local_ref_seq],
                           command="ln -s -f " + ref_seq + " " + local_ref_seq)
            job = concat_jobs([main_job, link_job], name="bismark_prepare_genome")
        else:  # In the off-chance that someone put the file in the output directory.
            job = main_job

        return [job]

    def bismark_align(self):
        """
        This step aligns trimmed reads to a bisulfite converted reference genome using Bismark. This create
        BAM files and will only be compressed if the input is also compressed (which usually will be the case).

        This step requires bismark_prepare_genome and the relevant trim_galore step.

        Input: Trimmed version of input files. (trimmed/*)
        Output: A BAM/SAM file in aligned/<sample_name>/*

        - [Set uneccessary files as removeable!]
        - [Get bismark into mugqic_pipelines modules]

        :return jobs: A list of jobs that needs to be executed in this step.
        :rtype list(Job):
        """

        # Multicore funtionality for Bismark is currently not supported when using the option --basename
        jobs = []
        for sample in self.samples:
            # Setup input/output directories and other configuration steps
            trim_prefix = os.path.join("trimmed", sample.name, sample.name)
            align_directory = os.path.join("aligned", sample.name)
            readset_sam = os.path.join(align_directory, sample.name + "_aligned_pe.sam.gz")
            run_type = sample.readsets[0].run_type

            # Specify input file name
            if run_type == "PAIRED_END":
                input_files = [trim_prefix + "_R1_val_1.fq.gz", trim_prefix + "_R2_val_2.fq.gz"]
            elif run_type == "SINGLE_END":
                input_files = [trim_prefix + "_R1_trimmed.fq.gz"]

            mkdir_job = Job(command="mkdir -p " + align_directory)
            job = concat_jobs([
                mkdir_job,
                Job(
                    input_files + ["Bisulfite_Genome"],
                    [readset_sam],
                    [["bismark_align", "module_bowtie2"],
                     ["bismark_align", "module_samtools"]],
                    command="""\
module load bismark/0.15
bismark -q --non_directional {other_options} --output_dir {directory} --basename {basename} . -1 {fastq1} -2 {fastq2}
    """.format(
                        directory=align_directory,
                        other_options=config.param("bismark_align", "other_options"),
                        fastq1=input_files[0],
                        fastq2=input_files[1] if run_type == "PAIRED_END" else "",
                        basename=sample.name + '_aligned'
                    )
                )], name="bismark_align." + sample.name)

            jobs.append(job)

        return jobs

    def bismark_deduplicate(self):
        """
        Calls the de-duplication module from Bismark. The module operates on the output alignment files from
        Bismark, creating a new output file (SAM default, BAM possible) with the suffix *.deduplicated.bam.
        The output file appears in the same directory as the input file, but this method will move the output file
        to its own directory.
        
        :return jobs: A list of jobs that needs to be executed in this step.
        :rtype list(Job):
        """

        jobs = []
        for sample in self.samples:
            work_dir = os.path.join('dedup', sample.name)
            in_file = os.path.join('aligned', sample.name, sample.name + '_aligned_pe.bam')
            out_file = os.path.join(work_dir, sample.name + '_aligned_pe.deduplicated.bam')
            run_type = sample.readsets[0].run_type
            protocol = sample.readsets[0].library

            if protocol == 'RRBS':  # Deduplication is not recommended for RRBS datatypes. Keep what we have
                job = Job([in_file], [out_file],
                          command="ln -s -f " + in_file + " " + out_file,
                          name="bismark_deduplication." + sample.name)
            else:
                job = concat_jobs([
                    Job([in_file],
                        [['bismark_deduplicate', 'module_samtools']],
                        command="""
module load bismark/0.15
deduplicate_bismark {type} --bam {other} {input}""".format(
                            type='--paired' if run_type == 'PAIRED' else '--single', input=in_file,
                            other=config.param('bismark_deduplicate', 'other_options', required=False)),
                        report_files=[os.path.join(work_dir, sample.name +
                                                   '_aligned_pe.deduplication_report.txt')]),
                    Job(output_files=[out_file],
                        command='mv -fu ' +
                                os.path.join('aligned', sample.name, sample.name +
                                             '_aligned_pe.deduplicate.bam') +
                                ' ' +
                                out_file)],
                    name='bismark_deduplication.' + sample.name)
            jobs.append(job)
        return jobs

    def bismark_methylation_caller(self):

        """
        This step extracts the methylation call for every single cytosine analyzed from the Bismark result files.
        The following input files are accepted:
            1.	Bismark result files from previous alignment step
            2.	BAM files (unsorted) from readset file

        Input: Merged sample files (merged/)
        Output: Methylation calls in BedGraph format. (methyl_calls/)

        :return jobs: A list of jobs that needs to be executed in this step.
        :rtype list(Job):
        """

        aligned_samples = []
        jobs = []

        for sample in self.samples:
            # Either select aligned sample from previous alignment step or aligned BAM/SAM files in readset file
            aligned_sample = os.path.join("aligned", sample.name, sample.name + "_aligned_pe.sam.gz")
            aligned_sample = self.select_input_files([[readset.bam for readset in sample.readsets], [
                os.path.join("aligned", sample.name, sample.name + "_aligned_pe.sam.gz")]])
            run_type = sample.readsets[0].run_type
            job = Job(
                aligned_sample,
                [os.path.join("methyl_calls", sample.name + "_aligned_pe.sam.bismark.cov.gz")],
                command="""\
        mkdir -p {directory}
        module load bismark/0.15
        bismark_methylation_extractor {library_type} {other_options} --output {directory} --multicore {cores} --bedGraph {sample}
        """.format(
                    directory="methyl_calls",
                    library_type="--paired-end" if run_type == "PAIRED_END" else "--single-end",
                    other_options=config.param("bismark_methylation_caller", "other_options"),
                    sample=" ".join(aligned_sample),
                    cores=config.param("bismark_methylation_caller", "cores", type="int") / 3
                ),
                name="bismark_methylation_caller." + sample.name)

            jobs.append(job)

        return jobs

    def differential_methylated_pos(self):

        """
        This step finds a list of differentially methylated CpG sites with respect to a categorical
        phenotype (controls vs. cases). The BedGraph files from the previous methylation calling step are first combined
        to a BSRaw object with the R package BiSeq. Then, the dmpFinder function from the R package minfi is used to
        compute a F-test statistic on the beta values for the assayed CpGs in each sample. A p-value is then returned
        for each site with the option of correcting them for multiple testing. Differential analysis is done for each
        contrast specified in the design file

        The values from the design files dictates how the samples are treated and compared.

        Input: Methylation data (methyl_calls/)
        Output: A CSV file in differential_methylated_positions/

        :return jobs: A list of jobs that needs to be executed in this step.
        :rtype list(Job):
        """

        jobs = []

        for contrast in self.contrasts:
            # Determine the control and case samples to include in the analysis from the contrast
            contrast_samples = [sample for sample in contrast.controls + contrast.treatments]
            cov_files = [os.path.join("methyl_calls", sample.name + "_aligned_pe.sam.bismark.cov.gz") for sample in
                         contrast_samples]
            sample_group = ["control" if sample in contrast.controls else "case" for sample in contrast_samples]
            dmps_file = os.path.join("differential_methylated_positions",
                                     contrast.name + "_RRBS_differential_methylated_pos.csv")

            job = Job(
                cov_files,
                [dmps_file],
                [
                    ["differential_methylated_pos", "module_R"],
                    ["differential_methylated_pos", "module_mugqic_R_packages"]
                ],
                command="""\
        mkdir -p {directory}
        R --no-save --no-restore <<-EOF
        suppressPackageStartupMessages(library(BiSeq))
        suppressPackageStartupMessages(library(minfi))
        rrbs <- readBismark(c{samples}, colData=DataFrame(group=factor(c{group}), row.names=c{sample_names}))
        rrbs.filtered <- rrbs[apply(totalReads(rrbs), 1, function(x) any(x > {coverage})),]
        beta <- methLevel(rawToRel(rrbs.filtered))

        #Use M values to do statistical tests because they are more reliable
        #dmpFinder does not work with M values that are 0 or INF so the beta values must be shifted slightly
        #Although there is no such thing as a beta value > 1, it will not matter in this step because only
        #the average beta values are shown to the user
        beta[beta == 0] = 0.0001
        beta[beta == 1] = 1.0001
        M <- log2(beta/(1-beta))

        dmp <- dmpFinder(M, pheno=colData(rrbs.filtered)[,"group"], type="categorical")
        dmp["pval"] <- p.adjust(dmp[,"pval"], method = "{padjust_method}")
        dmp <- dmp[dmp["pval"] < {pvalue},]["pval"]

        controls <- c{controls}
        cases <- c{cases}
        result = as.data.frame(rowRanges(rrbs.filtered))[1:4]
        result["Avg Control Beta"] = rowMeans(beta[,controls])
        result["Avg Case Beta"] = rowMeans(beta[,cases])
        result["Avg Delta Beta"] = result[,"Avg Case Beta"] - result[,"Avg Control Beta"]
        result <- merge(result, dmp, by=0)
        result <- result[abs(result["Avg Delta Beta"]) > {delta_beta_threshold}]

        write.csv(result, file="{dmps_file}", quote=FALSE, row.names=FALSE)

        EOF""".format(
                    directory=os.path.dirname(dmps_file),
                    samples=tuple(cov_files),
                    group=tuple(sample_group),
                    sample_names=tuple([sample.name for sample in contrast_samples]),
                    coverage=config.param("differential_methylated_pos", "read_coverage"),
                    controls=tuple([sample.name for sample in contrast.controls]),
                    cases=tuple([sample.name for sample in contrast.treatments]),
                    padjust_method=config.param("differential_methylated_pos", "padjust_method"),
                    pvalue=config.param("differential_methylated_pos", "pvalue", type="float"),
                    delta_beta_threshold=config.param("differential_methylated_pos", "delta_beta_threshold",
                                                      type="float"),
                    dmps_file=dmps_file
                ),
                name="differential_methylated_pos." + contrast.name)

            jobs.append(job)

        return jobs

    def differential_methylated_regions(self):

        """
        Similar to differential_methylated_positions, this step looks at methylation patterns on a larger, regional
        level. This step compares large-scale differences in methylation as opposed to comparing local methylation
        sites.

        Input: Methylation data (methyl_calls/)
        Output: A CSV file in differential_methylated_regions/

        :return jobs: A list of jobs that needs to be executed in this step.
        :rtype list(Job):
        """
        jobs = []

        for contrast in self.contrasts:
            # Determine the control and case samples to include in the analysis from the contrast
            contrast_samples = [sample for sample in contrast.controls + contrast.treatments]
            cov_files = [os.path.join("methyl_calls", sample.name + "_aligned_pe.sam.bismark.cov.gz") for sample in
                         contrast_samples]
            sample_group = ["control" if sample in contrast.controls else "case" for sample in contrast_samples]
            dmrs_file = os.path.join("differential_methylated_regions",
                                     contrast.name + "_RRBS_differential_methylated_regions.csv")

            job = Job(
                cov_files,
                [dmrs_file],
                [
                    ["differential_methylated_regions", "module_R"],
                    ["differential_methylated_regions", "module_mugqic_R_packages"]
                ],
                command="""\
        mkdir -p {directory}
        R --no-save --no-restore <<-EOF
        suppressPackageStartupMessages(library(bumphunter))
        suppressPackageStartupMessages(library(BiSeq))
        library(doParallel)
        registerDoParallel(cores={cores})

        rrbs <- readBismark(c{samples}, colData=DataFrame(group=c{group}, row.names=c{sample_names}))
        rrbs.filtered <- rrbs[apply(totalReads(rrbs), 1, function(x) any(x > {coverage})),]
        beta <- methLevel(rawToRel(rrbs.filtered))
        chr <- as.character(seqnames(rowRanges(rrbs.filtered)))
        pos <- start(ranges(rowRanges(rrbs.filtered)))
        pheno <- colData(rrbs.filtered)[,"group"]
        designM <- model.matrix(~pheno)

        dmrs <- bumphunterEngine(beta,
                                 chr=chr,
                                 pos=pos,
                                 design=designM,
                                 cutoff={delta_beta_threshold},
                                 pickCutoffQ=0.99,
                                 null_method=c("permutation","bootstrap"),
                                 smooth=FALSE,
                                 smoothFunction=locfitByCluster,
                                 B={permutations},
                                 verbose=TRUE,
                                 maxGap=500)

        dmrs <- na.omit(dmrs)

        write.csv(dmrs\$table, "{dmrs_file}", quote=FALSE, row.names=FALSE)

        EOF""".format(
                    directory=os.path.dirname(dmrs_file),
                    samples=tuple(cov_files),
                    group=tuple(sample_group),
                    coverage=config.param("differential_methylated_regions", "read_coverage", type="int"),
                    sample_names=tuple([sample.name for sample in contrast_samples]),
                    cores=config.param("differential_methylated_regions", "cores", type="int"),
                    delta_beta_threshold=config.param("differential_methylated_regions", "delta_beta_threshold",
                                                      type="float"),
                    permutations=config.param("differential_methylated_regions", "permutations", type="int"),
                    dmrs_file=dmrs_file
                ),
                name="differential_methylated_regions." + contrast.name)

            jobs.append(job)

        return jobs

    @property
    def steps(self):
        return [
            self.merge_fastq,
            self.trim_galore,
            self.bismark_prepare_genome,
            self.bismark_align,
            self.bismark_deduplicate,
            self.bismark_methylation_caller,
            self.differential_methylated_pos,
            self.differential_methylated_regions
        ]


if __name__ == '__main__':
    Episeq()
