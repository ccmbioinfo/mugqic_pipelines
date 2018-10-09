#!/usr/bin/env python

################################################################################
# Copyright (C) 2014, 2015 GenAP, McGill University and Genome Quebec Innovation Centre
#
# This file is part of MUGQIC Pipelines.
#
# MUGQIC Pipelines is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# MUGQIC Pipelines is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with MUGQIC Pipelines.  If not, see <http://www.gnu.org/licenses/>.
################################################################################

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
from bfx.design import *
from bfx.readset import *

from bfx import picard
from pipelines import common

from bfx import samtools_1_1
from bfx import defuse
from bfx import fusionmap
from bfx import tophat2
from bfx import integrate
from bfx import ericscript
from bfx import gunzip
from bfx import merge_fastq
from bfx import cff_conversion
from bfx import check_dna_support_before_next_exon
from bfx import merge_and_reannotate_cff_fusion
from bfx import repeat_filter
from bfx import fusion_stats
from bfx import validate_fusions
from bfx import delete_fastqs

import utils

from bfx import bedtools
from bfx import cufflinks
from bfx import differential_expression
from bfx import gq_seq_utils
from bfx import htseq
from bfx import metrics
from bfx import picard
from bfx import samtools
from bfx import star
from bfx import bvatools
from bfx import rmarkdown
from pipelines import common
import utils



log = logging.getLogger(__name__)

class RnaFusion(common.Illumina):
    """
    RNAFusion Pipeline
    ================
    The Gene Fusion pipeline identifies gene fusion events using RNA-seq FASTQ files.  
    
    Four separate tools detect fusion events: 
    [deFuse](https://sourceforge.net/p/defuse/wiki/DeFuse/), 
    [FusionMap](http://www.arrayserver.com/wiki/index.php?title=FusionMap), 
    [EricScript](https://sites.google.com/site/bioericscript/home), 
    and [INTEGRATE](https://sourceforge.net/p/integrate-fusion/wiki/Home/).
    
    Tophat2 is used to generate precursor files for the INTEGRATE fusion detection tool.
     
    The fusion detection results are combined into one common file (.cff) that gives information about gene fusions including gene names, 
    type of fusion (ex. read through vs. gene fusion), and the tools that identified each fusion event. 
    Additionally, if DNA sequencing information is available for the samples of interest, 
    the Gene Fusion Pipeline can check for DNA support of gene fusions detected from RNA. 
    
    The RNAseq pipeline requires a sampleinfo file to be provided, which is a tab-delimited document with each sample 
    as a line with information about sample disease. The first column gives sample name, second column gives the disease name,
    third column tells whether the sample comes from a tumor (TP) or normal (NT) tissue. Additional columns can give further 
    information about samples.

    In addition, a dnabam file must be provided, which gives the name of .bam file(s) associated with RNA-seq sample.
    If there is no DNA sequencing associated with the sample, provide the name of an empty file    

    For validation, an optional file containing fusion gene pairs can be provided with the --valfile flag. Used
    only if the validate_fusions step is being run, and assumes that the input sequence data contains the fusions provided in 
    the validation file. This step tests the effectiveness of the pipeline in detecting known fusions.

    """

    def __init__(self):
        # Add pipeline specific arguments
        self.argparser.add_argument("--sampleinfo", help="sample info file", type=file)
        self.argparser.add_argument("--dnabam", help="DNA bam list", type=file)
        # add optional fusion validation file for pipeline validation mode
        self.argparser.add_argument("--valfile", required=False, help="fusion validation set file", type=file)
        super(RnaFusion, self).__init__()


    def picard_sam_to_fastq(self):
        """
        Convert SAM/BAM files from the input readset file into FASTQ format
        if FASTQ files are not already specified in the readset file. Do nothing otherwise.
        rerwritten from common.Illumina.picard_sam_to_fastq, make directory for this step under result folder in case the orginal bam file directory is not writtable
        """
        jobs = []
        for readset in self.readsets:
            # If readset FASTQ files are available, skip this step
            if not readset.fastq1:
                if readset.cram:
                    # convert cram to bam then to fastq. fastq and bam are saved on localhd
                    out_bam = os.path.join("$TMPDIR", os.path.basename(readset.cram)+".bam")
                    cram2bam_job = samtools_1_1.view(readset.cram, out_bam)
                    if readset.run_type == "PAIRED_END":
                        out_dir = os.path.join("fusions", "picard_sam_to_fastq", readset.sample.name)
                        fastq1 = os.path.join(out_dir, os.path.basename(re.sub("\.bam$", ".pair1.fastq.gz", out_bam)))
                        fastq2 = os.path.join(out_dir, os.path.basename(re.sub("\.bam$", ".pair2.fastq.gz", out_bam)))
                    else:
                        raise Exception("Error: run type \"" + readset.run_type +
                        "\" is invalid for readset \"" + readset.name + "\" (should be PAIRED_END or SINGLE_END)!")

                    picard_job = picard.sam_to_fastq(out_bam, fastq1, fastq2)
                    job = concat_jobs([
                        Job(command="mkdir -p " + out_dir),
                        cram2bam_job,
                        picard_job
                    ], name= "picard_sam_to_fastq." + readset.name)
                    jobs.append(job)
                elif readset.bam:
                    if readset.run_type == "PAIRED_END":
                        out_dir = os.path.join("fusions", "picard_sam_to_fastq", readset.sample.name)
                        fastq1 = os.path.join(out_dir, os.path.basename(re.sub("\.bam$", ".pair1.fastq.gz", readset.bam)))
                        fastq2 = os.path.join(out_dir, os.path.basename(re.sub("\.bam$", ".pair2.fastq.gz", readset.bam)))
                    else:
                        raise Exception("Error: run type \"" + readset.run_type +
                        "\" is invalid for readset \"" + readset.name + "\" (should be PAIRED_END or SINGLE_END)!")

                    picard_job = picard.sam_to_fastq(readset.bam, fastq1, fastq2)
                    job = concat_jobs([
                        Job(command="mkdir -p " + out_dir),
                        picard_job
                    ], name= "picard_sam_to_fastq." + readset.name)
                    jobs.append(job)
                else:
                    raise Exception("Error: BAM file not available for readset \"" + readset.name + "\"!")
        return jobs
    def star(self):
        """
        The filtered reads are aligned to a reference genome. The alignment is done per readset of sequencing
        using the [STAR](https://code.google.com/p/rna-star/) software. It generates a Binary Alignment Map file (.bam).

        This step takes as input files:

        1. Trimmed FASTQ files if available
        2. Else, FASTQ files from the readset file if available
        3. Else, FASTQ output files from previous picard_sam_to_fastq conversion of BAM files
        """

        jobs = []
        project_index_directory = "reference.Merged"
        project_junction_file = os.path.join("alignment_1stPass", "AllSamples.SJ.out.tab")
        individual_junction_list=[]
        ######
        #pass 1 -alignment
        for readset in self.readsets:
            trim_file_prefix = os.path.join("trim", readset.sample.name, readset.name + ".trim.")
            alignment_1stPass_directory = os.path.join("alignment_1stPass", readset.sample.name, readset.name)
            individual_junction_list.append(os.path.join(alignment_1stPass_directory,"SJ.out.tab"))

            if readset.run_type == "PAIRED_END":
                candidate_input_files = [[trim_file_prefix + "pair1.fastq.gz", trim_file_prefix + "pair2.fastq.gz"]]
                if readset.fastq1 and readset.fastq2:
                    candidate_input_files.append([readset.fastq1, readset.fastq2])
                if readset.bam:
                    candidate_input_files.append([re.sub("\.bam$", ".pair1.fastq.gz", readset.bam), re.sub("\.bam$", ".pair2.fastq.gz", readset.bam)])
                [fastq1, fastq2] = self.select_input_files(candidate_input_files)
            elif readset.run_type == "SINGLE_END":
                candidate_input_files = [[trim_file_prefix + "single.fastq.gz"]]
                if readset.fastq1:
                    candidate_input_files.append([readset.fastq1])
                if readset.bam:
                    candidate_input_files.append([re.sub("\.bam$", ".single.fastq.gz", readset.bam)])
                [fastq1] = self.select_input_files(candidate_input_files)
                fastq2 = None
            else:
                raise Exception("Error: run type \"" + readset.run_type +
                "\" is invalid for readset \"" + readset.name + "\" (should be PAIRED_END or SINGLE_END)!")

            rg_platform = config.param('star_align', 'platform', required=False)
            rg_center = config.param('star_align', 'sequencing_center', required=False)

            job = star.align(
                reads1=fastq1,
                reads2=fastq2,
                output_directory=alignment_1stPass_directory,
                genome_index_folder=None,
                rg_id=readset.name,
                rg_sample=readset.sample.name,
                rg_library=readset.library if readset.library else "",
                rg_platform_unit=readset.run + "_" + readset.lane if readset.run and readset.lane else "",
                rg_platform=rg_platform if rg_platform else "",
                rg_center=rg_center if rg_center else ""
            )
            job.name = "star_align.1." + readset.name
            jobs.append(job)
        
        ######
        jobs.append(concat_jobs([
        #pass 1 - contatenate junction
        star.concatenate_junction(
            input_junction_files_list=individual_junction_list,
            output_junction_file=project_junction_file
        ),
        #pass 1 - genome indexing
        star.index(
            genome_index_folder=project_index_directory,
            junction_file=project_junction_file
        )], name = "star_index.AllSamples"))

        ######
        #Pass 2 - alignment
        for readset in self.readsets:
            trim_file_prefix = os.path.join("trim", readset.sample.name, readset.name + ".trim.")
            alignment_2ndPass_directory = os.path.join("alignment", readset.sample.name, readset.name)

            if readset.run_type == "PAIRED_END":
                candidate_input_files = [[trim_file_prefix + "pair1.fastq.gz", trim_file_prefix + "pair2.fastq.gz"]]
                if readset.fastq1 and readset.fastq2:
                    candidate_input_files.append([readset.fastq1, readset.fastq2])
                if readset.bam:
                    candidate_input_files.append([re.sub("\.bam$", ".pair1.fastq.gz", readset.bam), re.sub("\.bam$", ".pair2.fastq.gz", readset.bam)])
                [fastq1, fastq2] = self.select_input_files(candidate_input_files)
            elif readset.run_type == "SINGLE_END":
                candidate_input_files = [[trim_file_prefix + "single.fastq.gz"]]
                if readset.fastq1:
                    candidate_input_files.append([readset.fastq1])
                if readset.bam:
                    candidate_input_files.append([re.sub("\.bam$", ".single.fastq.gz", readset.bam)])
                [fastq1] = self.select_input_files(candidate_input_files)
                fastq2 = None
            else:
                raise Exception("Error: run type \"" + readset.run_type +
                "\" is invalid for readset \"" + readset.name + "\" (should be PAIRED_END or SINGLE_END)!")

            rg_platform = config.param('star_align', 'platform', required=False)
            rg_center = config.param('star_align', 'sequencing_center', required=False)

            job = star.align(
                reads1=fastq1,
                reads2=fastq2,
                output_directory=alignment_2ndPass_directory,
                genome_index_folder=project_index_directory,
                rg_id=readset.name,
                rg_sample=readset.sample.name,
                rg_library=readset.library if readset.library else "",
                rg_platform_unit=readset.run + "_" + readset.lane if readset.run and readset.lane else "",
                rg_platform=rg_platform if rg_platform else "",
                rg_center=rg_center if rg_center else "",
                create_wiggle_track=True,
                search_chimeres=True,
                cuff_follow=True,
                sort_bam=True
            )
            job.input_files.append(os.path.join(project_index_directory, "SAindex"))
 
            # If this readset is unique for this sample, further BAM merging is not necessary.
            # Thus, create a sample BAM symlink to the readset BAM.
            # remove older symlink before otherwise it raise an error if the link already exist (in case of redo)
            if len(readset.sample.readsets) == 1:
                readset_bam = os.path.join(alignment_2ndPass_directory, "Aligned.sortedByCoord.out.bam")
                sample_bam = os.path.join("alignment", readset.sample.name ,readset.sample.name + ".sorted.bam")
                job = concat_jobs([
                    job,
                    Job([readset_bam], [sample_bam], command="ln -s -f " + os.path.relpath(readset_bam, os.path.dirname(sample_bam)) + " " + sample_bam, removable_files=[sample_bam])])

            job.name = "star_align.2." + readset.name
            jobs.append(job)

        report_file = os.path.join("report", "RnaSeq.star.md")
        jobs.append(
            Job(
                [os.path.join("alignment", readset.sample.name, readset.name, "Aligned.sortedByCoord.out.bam") for readset in self.readsets],
                [report_file],
                [['star', 'module_pandoc']],
                command="""\
mkdir -p report && \\
pandoc --to=markdown \\
  --template {report_template_dir}/{basename_report_file} \\
  --variable scientific_name="{scientific_name}" \\
  --variable assembly="{assembly}" \\
  {report_template_dir}/{basename_report_file} \\
  > {report_file}""".format(
                    scientific_name=config.param('star', 'scientific_name'),
                    assembly=config.param('star', 'assembly'),
                    report_template_dir=self.report_template_dir,
                    basename_report_file=os.path.basename(report_file),
                    report_file=report_file
                ),
                report_files=[report_file],
                name="star_report")
        )
        return jobs
    def picard_merge_sam_files(self):
        """
        BAM readset files are merged into one file per sample. Merge is done using [Picard](http://broadinstitute.github.io/picard/).
        """

        jobs = []
        for sample in self.samples:
            # Skip samples with one readset only, since symlink has been created at align step
            if len(sample.readsets) > 1:
                alignment_directory = os.path.join("alignment", sample.name)
                inputs = [os.path.join(alignment_directory, readset.name, "Aligned.sortedByCoord.out.bam") for readset in sample.readsets]
                output = os.path.join(alignment_directory, sample.name + ".sorted.bam")

                job = picard.merge_sam_files(inputs, output)
                job.name = "picard_merge_sam_files." + sample.name
                jobs.append(job)
        return jobs

    def picard_sort_sam(self):
        """
        The alignment file is reordered (QueryName) using [Picard](http://broadinstitute.github.io/picard/). The QueryName-sorted bam files will be used to determine raw read counts.
        """

        jobs = []
        for sample in self.samples:
            alignment_file_prefix = os.path.join("alignment", sample.name, sample.name)

            job = picard.sort_sam(
                alignment_file_prefix + ".sorted.bam",
                alignment_file_prefix + ".QueryNameSorted.bam",
                "queryname"
            )
            job.name = "picard_sort_sam." + sample.name
            jobs.append(job)
        return jobs


    def gunzip_fastq(self):
        """
        Gunzip .fastq.gz files 
        """
        jobs = []
        for readset in self.readsets:
            trim_file_prefix = os.path.join("trim", readset.sample.name, readset.name + ".trim.")
            out_dir = os.path.join("fusions", "gunzip_fastq", readset.sample.name)
            # Find input readset FASTQs first from previous trimmomatic job, then from original FASTQs in the readset sheet
            if readset.run_type == "PAIRED_END":
                candidate_input_files = []
                if readset.fastq1 and readset.fastq2:
                    candidate_input_files.append([readset.fastq1, readset.fastq2])
                if readset.bam:
                    picard_dir = os.path.join("fusions", "picard_sam_to_fastq", readset.sample.name)
                    candidate_input_files.append([os.path.join(picard_dir, os.path.basename(re.sub("\.bam$", ".pair1.fastq.gz", readset.bam))), os.path.join(picard_dir, os.path.basename(re.sub("\.bam$", ".pair2.fastq.gz", readset.bam)))])
                if readset.cram:
                    picard_dir = os.path.join("fusions", "picard_sam_to_fastq", readset.sample.name)
                    candidate_input_files.append([os.path.join(picard_dir, os.path.basename(readset.cram)+".pair1.fastq.gz"), os.path.join(picard_dir, os.path.basename(readset.cram)+".pair2.fastq.gz")])
                [fastq1, fastq2] = self.select_input_files(candidate_input_files)
            else:
                raise Exception("Error: run type \"" + readset.run_type +
                "\" is invalid for readset \"" + readset.name + "\" (should be PAIRED_END)!")
            gunzip1_job = gunzip.gunzip_fastq(fastq1, out_dir)
            gunzip2_job = gunzip.gunzip_fastq(fastq2, out_dir)
            job = concat_jobs([
                Job(command="mkdir -p " + out_dir),
                gunzip1_job,
                gunzip2_job
            ], name="gunzip_fastq." + readset.sample.name + "." + readset.name)

            jobs.append(job)

        return jobs


    def merge_fastq(self):
        """
        Merge paired end fastqs of the same sample
        """
        jobs = []
        for sample in self.samples:
            if len(sample.readsets) > 1:
                input_dir = os.path.join("fusions", "gunzip_fastq", sample.name)
                fastq1_list = []
                fastq2_list = []
                for readset in sample.readsets:
                    if readset.bam:
                        fastq1 = os.path.join(input_dir, os.path.basename(re.sub("\.bam$", ".pair1.fastq", readset.bam)))
                        fastq2 = os.path.join(input_dir, os.path.basename(re.sub("\.bam$", ".pair2.fastq", readset.bam)))
                    if readset.fastq1:
                        if readset.fastq1.endswith(".gz"):
                            # input files are gzipped fastqs
                            fastq1 = os.path.join(input_dir, os.path.basename(re.sub("\.gz$", "", readset.fastq1)))
                            fastq2 = os.path.join(input_dir, os.path.basename(re.sub("\.gz$", "", readset.fastq2)))
                        else:
                            # input files are fastqs
                            fastq1 = os.path.join(input_dir, os.path.basename(readset.fastq1))
                            fastq2 = os.path.join(input_dir, os.path.basename(readset.fastq2))
                    fastq1_list.append(fastq1)
                    fastq2_list.append(fastq2)
                merge_fastq_job = merge_fastq.merge_fastq(fastq1_list, fastq2_list, input_dir)
                job = concat_jobs([
                    merge_fastq_job,
                ], name="merge_fastq." + sample.name)
                jobs.append(job)

        return jobs

    def select_input_fastq(self, sample):
        """
        Select input fastqs for fusion callers according to readset.
        This function is called in the gene fusion caller functions.
        """
        input_dir = os.path.join("fusions", "gunzip_fastq", sample.name)
        if len(sample.readsets) > 1: # sample has more than 1 readset, use merged fastq
            fastq1 = os.path.join(input_dir, "merged.pair1.fastq")
            fastq2 = os.path.join(input_dir, "merged.pair2.fastq")
        else:
            #input files are bams
            readset = sample.readsets[0]
            if readset.bam:
                fastq1 = os.path.join(input_dir, os.path.basename(re.sub("\.bam$", ".pair1.fastq", readset.bam)))
                fastq2 = os.path.join(input_dir, os.path.basename(re.sub("\.bam$", ".pair2.fastq", readset.bam)))
            if readset.fastq1:
                if readset.fastq1.endswith(".gz"):
                    # input files are gzipped fastqs
                    fastq1 = os.path.join(input_dir, os.path.basename(re.sub("\.gz$", "", readset.fastq1)))
                    fastq2 = os.path.join(input_dir, os.path.basename(re.sub("\.gz$", "", readset.fastq2)))
                else:
                    # input files are fastqs
                    fastq1 = os.path.join(input_dir, os.path.basename(readset.fastq1))
                    fastq2 = os.path.join(input_dir, os.path.basename(readset.fastq2))
            if readset.cram:
                fastq1 = os.path.join(input_dir, os.path.basename(readset.cram)+".pair1.fastq")
                fastq2 = os.path.join(input_dir, os.path.basename(readset.cram)+".pair2.fastq")
        print >> sys.stderr, fastq1
        print >> sys.stderr, fastq2
        return fastq1, fastq2


    def defuse(self):
        """
        Run Defuse to call gene fusions
        """
        jobs = []
        for sample in self.samples:
            fastq1, fastq2 = self.select_input_fastq(sample)
            out_dir = os.path.join("fusions", "defuse", sample.name)
            defuse_job = defuse.defuse(fastq1, fastq2, out_dir)
            job = concat_jobs([
                Job(command="mkdir -p " + out_dir),
                defuse_job
            ], name="defuse." + sample.name)

            jobs.append(job)

        return jobs


    def fusionmap(self):
        """
        Run FusionMap to call gene fusions
        """
        jobs = []
        for sample in self.samples:
            fastq1, fastq2 = self.select_input_fastq(sample)
            out_dir = os.path.join("fusions", "fusionmap", sample.name)
            fusionmap_job = fusionmap.fusionmap(fastq1, fastq2, out_dir)
            job = concat_jobs([
                Job(command="mkdir -p " + out_dir),
                fusionmap_job,
                Job(command="ls " + out_dir + "/02_RNA*")
            
            ], name="fusionmap." + sample.name)

            jobs.append(job)

        return jobs


    def ericscript(self):
        """
        Run EricScript to call gene fusions
        """
        jobs = []
        for sample in self.samples:
            fastq1, fastq2 = self.select_input_fastq(sample)
            out_dir = os.path.join("fusions", "ericscript", sample.name)
            ericscript_job = ericscript.ericscript(fastq1, fastq2, out_dir)
            job = concat_jobs([
                Job(command="mkdir -p " + out_dir),
                Job(command="rm -r " + out_dir),
                ericscript_job
            ], name="ericscript." + sample.name)

            jobs.append(job)

        return jobs    
    
    def tophat2(self):
        """
        Run Tophat2 for Integrate. Determines accepted hits and unmapped reads, and outputs 
        corresponding .bam files required as input files for integrate step.
        """
        jobs = []
        for sample in self.samples:
            fastq1, fastq2 = self.select_input_fastq(sample)
            out_dir = os.path.join(self.output_dir, "fusions", "tophat2", sample.name)
            tophat2_job = tophat2.tophat2(fastq1, fastq2, out_dir)
            job = concat_jobs([
                Job(command="mkdir -p " + out_dir),
                tophat2_job
            ], name="tophat2." + sample.name)

            jobs.append(job)

        return jobs    

    def integrate(self):
        """
        Run Integrate to call gene fusions
        """
        jobs = []
        for sample in self.samples:
            input_dir = os.path.join("fusions", "tophat2", sample.name)
            accepted_bam = os.path.join(self.output_dir, input_dir, "accepted_hits.bam")
            unmapped_bam = os.path.join(self.output_dir, input_dir, "unmapped.bam")

            out_dir = os.path.join("fusions", "integrate", sample.name)
            integrate_job = integrate.integrate(accepted_bam, unmapped_bam, out_dir)
            job = concat_jobs([
                Job(command="mkdir -p " + out_dir),
                Job(command="cd " + out_dir),
                integrate_job,
                Job(command="cd -")
            ], name="integrate." + sample.name)

            jobs.append(job)

        return jobs
    
    def integrate_make_result_file(self):
        """
        Merge infomation from breakpoints.tsv and reads.txt
        """
        jobs = []
        for sample in self.samples:
            input_dir = os.path.join("fusions", "integrate", sample.name)

            make_result_job = integrate.make_result_file(input_dir)
            job = concat_jobs([
                make_result_job
            ], name="integrate_make_result." + sample.name)

            jobs.append(job)

        return jobs    
        

    def convert_fusion_results_to_cff(self):
        """
        Convert fusion results of all 4 gene fusion callers to cff format
        """
        jobs = []
        out_dir = os.path.join("fusions", "cff")
        job_list = [Job(command="mkdir -p " + out_dir)]
        sampleinfo_file = os.path.relpath(self.args.sampleinfo.name, self.output_dir)
        for sample in self.samples:
            defuse_result = os.path.join("fusions", "defuse", sample.name, "results.filtered.tsv")
            fusionmap_result = os.path.join("fusions", "fusionmap", sample.name, "02_RNA.FusionReport.txt")
            ericscript_result = os.path.join("fusions", "ericscript", sample.name, "fusion.results.filtered.tsv")
            integrate_result = os.path.join("fusions", "integrate", sample.name, "breakpoints.cov.tsv")

            tool_results = [("defuse", defuse_result), ("fusionmap", fusionmap_result), ("ericscript", ericscript_result), ("integrate", integrate_result)]
            """
            sample_type = ""
            for contrast in self.contrasts:
                if sample in contrast.controls:
                    sample_type = "Normal"
                elif sample in contrast.treatments:
                    sample_type = "Tumor"
                if sample_type:
                    disease_name = contrast.name
                    break    
            if not sample_type:
                raise Exception("Error: sample " + sample.name + " not found in design file " + self.args.design.name)
            """    
            for tool, result_file in tool_results:
                job = cff_conversion.cff_convert(sample.name, result_file, sampleinfo_file, tool, out_dir)
                job.command = job.command.strip()
                job_list.append(job)
        job = concat_jobs(job_list, name="cff_conversion")
        jobs.append(job)
        return jobs

    def merge_and_reannotate_cff_fusion(self):
        """
        Merge all cff files into one single file and reannotate it with given annotation files
        """
        jobs = []
        cff_files = []
        cff_dir = os.path.join("fusions", "cff")
        out_dir = os.path.join("fusions", "cff")
        tool_list = ["defuse", "fusionmap", "ericscript", "integrate"]
        for tool in tool_list:
            cff_files.extend([os.path.join(cff_dir, sample.name+"."+tool+".cff") for sample in self.samples])
        
        reann_job = merge_and_reannotate_cff_fusion.merge_and_reannotate_cff_fusion(cff_files, out_dir)
        
        job = concat_jobs([
            reann_job    
        ], name="merge_and_reannotate_cff_fusion")

        jobs.append(job)
        return jobs
        
    def check_dna_support_before_next_exon(self):
        """
        Check DNA support (pair clusters) until the start of next exon/utr
        """
        jobs = []
        dna_bam_list = os.path.abspath(self.args.dnabam.name)
        tmp_dir = os.path.join("fusions", "tmp")
        reann_file = os.path.join("fusions", "cff", "merged.cff.reann")
        dna_supp_job = check_dna_support_before_next_exon.check_dna_support_before_next_exon(reann_file, dna_bam_list, tmp_dir)
        job = concat_jobs([
            Job(command="mkdir -p " + tmp_dir),
            dna_supp_job    
        ], name="check_dna_support_before_next_exon")

        jobs.append(job)
        return jobs


    def repeat_filter(self):
        """
        Filter fusions with repetitive boundary sequences by realigning a certain length of sequnces with BWA
        """
        jobs = []
        out_dir = os.path.join("fusions", "cff")
        cff = os.path.join(out_dir, "merged.cff.reann.dnasupp")
        job = repeat_filter.repeat_filter(cff, out_dir)
        
        job = concat_jobs([
            job    
        ], name="repeat_filter")

        jobs.append(job)
        return jobs

    def cluster_reann_dnasupp_file(self):
        """
        Reannotate DNA support (pair clusters) file. This step generates the final category/cluster file,
        merged.cff.reann.dnasupp.bwafilter.30.cluster, which has the following columns:
            cluster_type, gene1, gene2, max_split_cnt, max_span_cnt, sample_type, disease, tools, inferred_fusion_type,
            gene1_on_bnd, gene1_close_to_bnd, gene2_on_bnd, gene2_close_to_bnd, dna_supp, samples
  
        Gene_Cluster    GTF2I   PGS1    2   -1  Tumor   LIS fusionmap   GeneFusion  True    True    True    True    -1  LIS_S8_L001

        """
        jobs = []
        out_dir = os.path.join("fusions", "cff")
        cluster_job = merge_and_reannotate_cff_fusion.cluster_reann_dnasupp_file(out_dir)
        
        job = concat_jobs([
            cluster_job    
        ], name="cluster_reann_dnasupp_file")

        jobs.append(job)
        return jobs

    
    def fusion_stats(self):
        """
        Genereates a file containing statistics about detected fusions.
        """
        jobs = []
        in_dir = os.path.join("fusions", "cff")
        out_dir= os.path.join("fusions", "fusion_stats")
        sampleinfo_file = os.path.relpath(self.args.sampleinfo.name, self.output_dir)

        fusion_stats_job = fusion_stats.fusion_stats(in_dir, out_dir, sampleinfo_file)
        
        job = concat_jobs([
            Job(command="mkdir -p " + out_dir),
            fusion_stats_job
        ], name="fusion_stats")

        jobs.append(job)
        return jobs

    def validate_fusions(self):
        """
        Compares the pipeline output in merged.cff.reann.dnasupp.bwafilter.30.cluster with the predetermined
        fusion gene test file. Outputs statistics about the detected gene fusions. This 
        step should be run only with a test .bam/.fastq file, in order to validate fusions that are known to be
        present in this sequence data.
        Requires that the --valfile flag is used, and that the input file has the 5' and 3' fusion genes in the first
        and second column, respectively.
        """
        #output_fusions = os.path.join("fusions", "cff", "merged.cff.reann.dnasupp.bwafilter.30.cluster")
        in_dir = os.path.join("fusions", "cff")
        out_dir = os.path.join("fusions", "validate_fusions")
        jobs = []
        
        validate_fusions_job = validate_fusions.validate_fusions(in_dir, out_dir, self.args.valfile.name)

        job = concat_jobs([
            Job(command="mkdir -p " + out_dir),
            validate_fusions_job
        ], name="validate_fusions")

        jobs.append(job)
        return jobs

    def delete_fastqs(self):
        """
        Delete fastqs when all callers' jobs are finished
        """
        jobs = []
        for sample in self.samples:
            defuse_result = os.path.join("fusions", "defuse", sample.name, "results.filtered.tsv")
            fusionmap_result = os.path.join("fusions", "fusionmap", sample.name, "02_RNA.FusionReport.txt")
            ericscript_result = os.path.join("fusions", "ericscript", sample.name, "fusion.results.filtered.tsv")
            integrate_result = os.path.join("fusions", "integrate", sample.name, "breakpoints.cov.tsv")
            result_file_list = [defuse_result, fusionmap_result, ericscript_result, integrate_result]
            del_job = delete_fastqs.delete_fastqs(sample.name, result_file_list)
            job = concat_jobs([
                Job(command="mkdir -p delete_fastqs"),
                del_job
            ], name="delete_fastqs." + sample.name)
            jobs.append(job)
        return jobs

    @property
    def steps(self):
        return [
            self.picard_sam_to_fastq,
            self.trimmomatic,
            self.merge_trimmomatic_stats,
            self.star,
            self.picard_merge_sam_files,
            self.picard_sort_sam,
            self.gunzip_fastq,
            self.merge_fastq,
            self.defuse,
            self.fusionmap,
            self.ericscript,
            self.tophat2,
            self.integrate,
            self.integrate_make_result_file,
            self.convert_fusion_results_to_cff,
            self.merge_and_reannotate_cff_fusion,
            self.check_dna_support_before_next_exon,
            self.repeat_filter,
            self.cluster_reann_dnasupp_file,
            self.fusion_stats,
            self.validate_fusions,
            self.delete_fastqs
        ]

if __name__ == '__main__':
    RnaFusion()
