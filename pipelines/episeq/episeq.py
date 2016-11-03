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

    def bismark_prepare_genome(self):
        """
        Bismark requires a processed reference genome to compare with the epigenome.
        """
        ref_seq = config.param('bismark_prepare_genome', 'genome_file', type='filepath')
        local_ref_seq = os.path.join('bismark_prepare_genome', os.path.basename(ref_seq))
        output_idx = "Bisulfite_Genome"

        if not os.path.isfile(local_ref_seq):
            run_job = Job([ref_seq], [output_idx],
                          module_entries=[["bismark_prepare_genome", "module_bowtie2"],
                                          ["bismark_prepare_genome", "module_samtools"]],
                          command="""\
    cp {src} .
    module load bismark/0.16.3
    bismark_genome_preparation {verbose} --bowtie2 .""".format(
                              src=ref_seq,
                              verbose='--verbose' if config.param('bismark_genome_preparation', 'verbose_logging',
                                                                  required=False, type='boolean') else ''),
                          name="bismark_prepare_genome")
        else:
            run_job = Job([ref_seq], [output_idx],
                          module_entries=[["bismark_prepare_genome", "module_bowtie2"],
                                          ["bismark_prepare_genome", "module_samtools"]],
                          command="""\
    module load bismark/0.16.3
    bismark_genome_preparation {verbose} --bowtie2 .""".format(
                              verbose='--verbose' if config.param('bismark_genome_preparation', 'verbose_logging',
                                                                  required=False, type='boolean') else ''),
                          name="bismark_prepare_genome")

        return [run_job]

    def trim_galore(self):
        """
        This step trims raw FASTQ files for quality control using Trim Galore!
        This is a pre-proccessing step to ensure quality control.

        :return jobs: A list of jobs that needs to be executed in this step.
        :rtype list(Job):
        """

        jobs = []
        for sample in self.samples:
            for readset in sample.readsets:
                if readset.bam and not readset.fastq1:
                    continue  # There's still a bam file to process later
                elif not (readset.bam or readset.fastq1 or readset.fastq2):
                    raise ValueError("There is no data files in this readset: " + readset.name + "!\n")
                # Get metadata and names
                run_type = readset.run_type
                protocol = readset.library
                trim_directory = os.path.join("trimmed", sample.name, readset.name)
                fq1_out = os.path.join(trim_directory, readset.name)
                fq2_out = os.path.join(trim_directory, readset.name)
                output_files = []

                # Trim Galoree has no built in option to change the filenames of the output
                # Below are the default output names when running in paired or single mode
                if run_type == "PAIRED_END":
                    input_files = [readset.fastq1, readset.fastq2]
                    output_files = [fq1_out + "_1_val_1.fq.gz", fq2_out + "_2_val_2.fq.gz"]
                elif run_type == "SINGLE_END":
                    input_files = [readset.fastq1]
                    output_files = [fq1_out + "_trimmed.fq.gz"]

                mkdir_job = Job(command="mkdir -p " + trim_directory)
                job = concat_jobs([
                    mkdir_job,
                    Job(
                        input_files,
                        output_files,
                        [['trim_galore', 'module_fastqc']],
                        command="""\
    module load trim_galore/0.4.1
    module load cutadapt/1.10
    trim_galore {protocol} {library_type} {non_directional} {other} --output_dir {directory} {fastq}
            """.format(
                            library_type="--paired" if run_type == "PAIRED_END" else "",
                            protocol='--rrbs' if protocol == 'RRBS' else '',
                            non_directional='--non_directional' if run_type == 'PAIRED_END' and protocol == 'RRBS'
                            else '',
                            other=config.param("trim_galore", "other_options"),
                            directory=trim_directory,
                            fastq=' '.join(input_files)
                        )
                    )], name="trim_galore." + readset.name)
                jobs.append(job)
        return jobs

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
        # TODO: Add case for bam file input.
        # Multicore funtionality for Bismark is currently not supported when using the option --basename
        jobs = []
        for sample in self.samples:
            for readset in sample.readsets:
                run_type = readset.run_type
                protocol = readset.library
                trim_prefix = os.path.join("trimmed", sample.name, readset.name)
                align_directory = os.path.join("aligned", sample.name)
                readset_base = os.path.join(align_directory, readset.name)

                if run_type == "PAIRED_END":
                    input_files = [os.path.join(trim_prefix, readset.name + "_1_val_1.fq.gz"),
                                   os.path.join(trim_prefix, readset.name + "_2_val_2.fq.gz")]
                    readset_sam = readset_base + "_aligned_pe.bam"
                elif run_type == "SINGLE_END":
                    input_files = [os.path.join(trim_prefix, readset.name + "_trimmed.fq.gz")]
                    readset_sam = readset_base + "_aligned.bam"

                # Case when only a bam file is given for a readset
                if not readset.fastq1:
                    if readset.bam:
                        continue
                    else:
                        raise IOError("""{readset} has no input files!""".format(readset=readset.name))
                if not readset_sam:
                    raise AttributeError("Unknown run_type: " + run_type + ". Unknown file output name.")

                mkdir_job = Job(command="mkdir -p " + align_directory)
                job = concat_jobs([
                    mkdir_job,
                    Job(
                        input_files + ["Bisulfite_Genome"],
                        [readset_sam],
                        [["bismark_align", "module_bowtie2"],
                         ["bismark_align", "module_samtools"]],
                        # ['bismark_align', 'module_perl']],
                        # Do not import. Bismark is expecting /usr/bin/perl. Will raise errors otherwise.
                        command="""\
        module load bismark/0.16.3
        bismark -q {directional} {other} --output_dir {directory} --basename {basename} --genome_folder . {input}
        """.format(
                            directory=align_directory,
                            other=config.param("bismark_align", "other_options"),
                            input="-1 {fastq1} -2 {fastq2}".format(fastq1=input_files[0], fastq2=input_files[1]) if
                            run_type == "PAIRED_END" else "--single_end {fastq1}".format(fastq1=input_files[0]),
                            directional='--non_directional' if protocol == 'RRBS' else '',
                            basename=readset.name + '_aligned'
                        )
                    )], name="bismark_align." + readset.name)
                jobs.append(job)
        return jobs

    def picard_merge_sam_files(self):
        """

        :return:
        :rtype:
        """

        jobs = []
        for sample in self.samples:
            # Bam files from pipeline (FASTQ)
            processed_fastq_pe = [os.path.join('aligned', sample.name, readset.name + "_aligned_pe.bam") for
                                  readset in sample.readsets if readset.run_type == 'PAIRED_END' and not readset.bam]
            processed_fastq_se = [os.path.join('aligned', sample.name, readset.name + "_aligned.bam") for
                                  readset in sample.readsets if readset.run_type == 'SINGLE_END' and not readset.bam]
            processed_fastq = processed_fastq_pe + processed_fastq_se
            # Bam files from user, if specified by readset file. Not exclusive with having fastq
            listed_bam_files = [readset.bam for readset in sample.readsets if readset.bam != '' and not readset.fastq1]

            input_files = processed_fastq + listed_bam_files  # All bam files that belong to the sample
            merge_prefix = 'merged'
            output_bam = os.path.join(merge_prefix, sample.name + '.merged.bam')

            mkdir_job = Job(command='mkdir -p ' + merge_prefix)

            # I want to use Picard Tools v2.0.1, which has a different syntax than v1.x
            if len(input_files) > 1:
                picard_v2 = Job(
                    input_files,
                    [output_bam, re.sub("\.([sb])am$", ".\\1ai", output_bam)],
                    [
                        ['picard_merge_sam_files', 'module_java'],
                        ['picard_merge_sam_files', 'module_picard']
                    ],
                    command="""\
    java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} \\
      -jar $PICARD_HOME/picard.jar MergeSamFiles \\
      VALIDATION_STRINGENCY=SILENT \\
      TMP_DIR={tmp_dir} \\
      {inputs} \\
      OUTPUT={output} \\
      USE_THREADING=true \\
      SORT_ORDER=queryname \\
      MAX_RECORDS_IN_RAM={max_records_in_ram}""".format(
                        tmp_dir=config.param('picard_merge_sam_files', 'tmp_dir'),
                        java_other_options=config.param('picard_merge_sam_files', 'java_other_options'),
                        ram=config.param('picard_merge_sam_files', 'ram'),
                        inputs=" \\\n  ".join(["INPUT=" + in_put for in_put in input_files]),
                        output=output_bam,
                        max_records_in_ram=config.param('picard_merge_sam_files', 'max_records_in_ram', type='int')),
                    removable_files=[output_bam, re.sub("\.([sb])am$", ".\\1ai", output_bam)],
                    local=config.param('picard_merge_sam_files', 'use_localhd', required=False))

                job = concat_jobs([mkdir_job, picard_v2], name="picard_merge_sam_files." + sample.name)

            elif len(input_files) == 1:
                readset_bam = input_files[0]
                if os.path.isabs(readset_bam):
                    target_readset_bam = readset_bam
                else:
                    target_readset_bam = os.path.relpath(readset_bam, merge_prefix)
                job = concat_jobs([
                    mkdir_job,
                    Job([readset_bam], [output_bam], command="ln -s -f " + target_readset_bam + " " + output_bam,
                        removable_files=[output_bam])],
                    name="symlink_readset_sample_bam." + sample.name)
            else:
                raise ValueError('Sample ' + sample.name + ' has no readsets!')
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

        jobs = []

        for sample in self.samples:
            # Either select aligned sample from previous alignment step or aligned BAM/SAM files in readset file
            merged_sample = self.select_input_files([[readset.bam for readset in sample.readsets], [
                os.path.join("merged", sample.name + ".merged.bam")]])
            run_type = sample.readsets[0].run_type
            job = Job(
                merged_sample,
                [os.path.join("methyl_calls", sample.name, sample.name + ".merged.bismark.cov.gz")],
                [['bismark_methylation_caller', 'module_samtools']],
                command="""\
        mkdir -p {directory}
        module load bismark/0.16.3
        bismark_methylation_extractor {library_type} {other} --output {directory} --bedGraph --genome_folder . {sample}
        """.format(
                    directory=os.path.join("methyl_calls",  sample.name),
                    library_type="--paired-end" if run_type == "PAIRED_END" else "--single-end",
                    other=config.param("bismark_methylation_caller", "other_options"),
                    sample=" ".join(merged_sample)),
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
            cov_files = [os.path.join("methyl_calls", sample.name, sample.name + ".merged.bismark.cov.gz") for
                         sample in contrast_samples]
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

EOF
                """.format(
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
            cov_files = [os.path.join("methyl_calls", sample.name, sample.name + ".merged.bismark.cov.gz") for sample in
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

EOF
                """.format(
                    directory=os.path.dirname(dmrs_file),
                    samples=tuple(cov_files),
                    group=tuple(sample_group),
                    coverage=config.param("differential_methylated_regions", "read_coverage", type="int"),
                    sample_names=tuple([sample.name for sample in contrast_samples]),
                    cores=config.param('bismark_methylation_caller', 'cluster_cpu').split('=')[-1],
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
        """

        :return: list
        """
        return [
            self.bismark_prepare_genome,
            self.trim_galore,
            self.bismark_align,
            self.picard_merge_sam_files,
            self.bismark_methylation_caller,
            self.differential_methylated_pos,
            self.differential_methylated_regions
        ]


if __name__ == '__main__':
    Episeq()
