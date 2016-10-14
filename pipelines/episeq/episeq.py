#!/usr/bin/env python

# Epigenetics pipeline for RRBS/WGSB data

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
    The Episeq pipeline takes FASTQ or BAM files (unsorted) as input
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
    module load bismark/0.15
    bismark_genome_preparation --verbose --bowtie2 .""".format(src=ref_seq),
                          name="bismark_prepare_genome")
        else:
            # Run bismark
            run_job = Job([ref_seq], [output_idx],
                          module_entries=[["bismark_prepare_genome", "module_bowtie2"],
                                          ["bismark_prepare_genome", "module_samtools"]],
                          command="""\
    module load bismark/0.15
    bismark_genome_preparation --verbose --bowtie2 .""",
                          name="bismark_prepare_genome")

        return [run_job]

    def trim_galore(self):

        """
        This step trims raw FASTQ files for quality control using Trim Galore! and runs fastqc at the end,
        courtesy of trim galore.
        """

        jobs = []
        for sample in self.samples:
            for readset in sample.readsets:
                run_type = readset.run_type
                protocol = readset.library
                trim_directory = os.path.join("trimmed", sample.name, readset.name)
                fq1_out = os.path.join(trim_directory, os.path.splitext(readset.fastq1))
                fq2_out = os.path.join(trim_directory, os.path.splitext(readset.fastq2))
                # trimmed/sample.name/readset.name/sample.name_readset.name_R1_val_1.fq.gz
                output_files = []

                # Trim Galoree has no built in option to change the filenames of the output
                # Below are the default output names when running in paired or single mode
                if run_type == "PAIRED_END":
                    input_files = [readset.fastq1, readset.fastq2]
                    output_files = [fq1_out + "_R1_val_1.fq.gz", fq2_out + "_R2_val_2.fq.gz"]
                elif run_type == "SINGLE_END":
                    input_files = [readset.fastq1]
                    output_files = [fq1_out + "_R1_trimmed.fq.gz"]

                mkdir_job = Job(command="mkdir -p " + trim_directory)
                job = concat_jobs([
                    mkdir_job,
                    Job(
                        input_files,
                        output_files,
                        ['trim_galore', 'module_fastqc'],
                        command="""\
    module load trim_galore/0.4.1
    module load cutadapt/1.10
    trim_galore {protocol} {library_type} {non_directional} {other_options} --output_dir {directory} --fastqc\
    --fastqc_args "--outdir {directory} --dir {directory} --threads" {fastq1} {fastq2}
        """.format(
                            library_type="--paired" if run_type == "PAIRED_END" else "",
                            protocol='--rrbs' if protocol == 'RRBS' else '',
                            non_directional='--non_directional' if run_type == 'PAIRED_END' and protocol == 'RRBS'
                            else '',
                            other_options=config.param("trim_galore", "other_options"),
                            directory=trim_directory,
                            fastq1=input_files[0],
                            fastq2=input_files[1] if run_type == "PAIRED_END" else ""
                        )
                    )], name="trim_galore." + readset.name)
                jobs.append(job)
        return jobs

    def bismark_align(self):

        """
        This step aligns trimmed reads to a bisulfite converted reference genome using Bismark [link]
        """

        # Multicore funtionality for Bismark is currently not supported when using the option --basename
        jobs = []
        for sample in self.samples:
            trim_prefix = os.path.join("trimmed", sample.name, sample.name)
            align_directory = os.path.join("aligned", sample.name)
            readset_sam = os.path.join(align_directory, sample.name + "_aligned_pe.sam.gz")
            run_type = sample.readsets[0].run_type

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

    def bismark_methylation_caller(self):

        """
        This step extracts the methylation call for every single cytosine analyzed from the Bismark result files.
        The following input files are accepted:
            1.	Bismark result files from previous alignment step
            2.	BAM files (unsorted) from readset file
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
            self.bismark_methylation_caller,
            self.differential_methylated_pos,
            self.differential_methylated_regions
        ]


if __name__ == '__main__':
    Episeq()
