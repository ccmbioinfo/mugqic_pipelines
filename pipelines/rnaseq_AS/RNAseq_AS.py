#!/usr/bin/env python

# This is the RNAseq pipeline with modifications to allow alternative splicing
# from tools such as Miso, Vast Tools, Cufflinks and JunctionSeq.

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
from bfx import miso_funcs
from bfx import vtools
from bfx import jctseq
from pipelines import common
import utils

log = logging.getLogger(__name__)

class RnaSeq(common.Illumina):
    """
    RNA-Seq Pipeline
    ================

    The standard MUGQIC RNA-Seq pipeline is based on the use of the [STAR aligner](https://code.google.com/p/rna-star/)
    to align reads to the reference genome. These alignments are used during
    downstream analysis to determine genes and transcripts differential expression. The
    [Cufflinks](http://cufflinks.cbcb.umd.edu/) suite is used for the transcript analysis whereas
    [DESeq](http://bioconductor.org/packages/release/bioc/html/DESeq.html) and
    [edgeR](http://bioconductor.org/packages/release/bioc/html/edgeR.html) are used for the gene analysis.

    The RNAseq pipeline requires to provide a design file which will be used to define group comparison
    in the differential analyses. The design file format is described
    [here](https://bitbucket.org/mugqic/mugqic_pipelines/src#markdown-header-design-file)

    The differential gene analysis is followed by a Gene Ontology (GO) enrichment analysis.
    This analysis use the [goseq approach](http://bioconductor.org/packages/release/bioc/html/goseq.html).
    The goseq is based on the use of non-native GO terms (see details in the section 5 of
    [the corresponding vignette](http://bioconductor.org/packages/release/bioc/vignettes/goseq/inst/doc/goseq.pdf).

    Finally, a summary html report is automatically generated by the pipeline at the end of the analysis.
    This report contains description
    of the sequencing experiment as well as a detailed presentation of the pipeline steps and results.
    Various Quality Control (QC) summary statistics are included in the report and additional QC analysis
    is accessible for download directly through the report. The report includes also the main references
    of the software tools and methods used during the analysis, together with the full list of parameters
    that have been passed to the pipeline main script.

    An example of the RNA-Seq report for an analysis on Public Corriel CEPH B-cell is available for illustration
    purpose only: [RNA-Seq report](http://gqinnovationcenter.com/services/bioinformatics/tools/rnaReport/index.html).

    [Here](https://bitbucket.org/mugqic/mugqic_pipelines/downloads/MUGQIC_Bioinfo_RNA-Seq.pptx) is more
    information about the RNA-Seq pipeline that you may find interesting.
    """

    def __init__(self):
        # Add pipeline specific arguments
        self.argparser.add_argument("-d", "--design", help="design file", type=file)
        super(RnaSeq, self).__init__()

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

    def picard_mark_duplicates(self):
        """
        Mark duplicates. Aligned reads per sample are duplicates if they have the same 5' alignment positions
        (for both mates in the case of paired-end reads). All but the best pair (based on alignment score)
        will be marked as a duplicate in the BAM file. Marking duplicates is done using [Picard](http://broadinstitute.github.io/picard/).
        """

        jobs = []
        for sample in self.samples:
            alignment_file_prefix = os.path.join("alignment", sample.name, sample.name + ".sorted.")

            job = picard.mark_duplicates(
                [alignment_file_prefix + "bam"],
                alignment_file_prefix + "mdup.bam",
                alignment_file_prefix + "mdup.metrics"
            )
            job.name = "picard_mark_duplicates." + sample.name
            jobs.append(job)
        return jobs

    def bam_hard_clip(self):
        """
        Generate a hardclipped version of the bam for the toxedo suite which doesn't support this official sam feature.
        """
        
        jobs = []
        for sample in self.samples:
            alignment_input = os.path.join("alignment", sample.name, sample.name + ".sorted.mdup.bam")
            alignment_output = os.path.join("alignment", sample.name, sample.name + ".sorted.mdup.hardClip.bam")
            job=pipe_jobs([
                samtools.view(
                    alignment_input,
                    None,
                    "-h"
                ),
                Job(
                    [None],
                    [alignment_output],
                    # awk to transform soft clip into hard clip for tuxedo suite
                    command="""\
awk 'BEGIN {{OFS="\\t"}} {{if (substr($1,1,1)=="@") {{print;next}}; split($6,C,/[0-9]*/); split($6,L,/[SMDIN]/); if (C[2]=="S") {{$10=substr($10,L[1]+1); $11=substr($11,L[1]+1)}}; if (C[length(C)]=="S") {{L1=length($10)-L[length(L)-1]; $10=substr($10,1,L1); $11=substr($11,1,L1); }}; gsub(/[0-9]*S/,"",$6); print}}' """.format()
                ),
                samtools.view(
                    "-",
                    alignment_output,
                    "-hbS"
                ),
            ])
            job.name="tuxedo_hard_clip."+ sample.name
            jobs.append(job)
        return jobs
    

    def picard_rna_metrics(self):
        """
        Computes a series of quality control metrics using both CollectRnaSeqMetrics and CollectAlignmentSummaryMetrics functions
        metrics are collected using [Picard](http://broadinstitute.github.io/picard/).
        """
        
        jobs = []
        reference_file = config.param('picard_rna_metrics', 'genome_fasta', type='filepath')
        for sample in self.samples:
                alignment_file = os.path.join("alignment", sample.name, sample.name + ".sorted.mdup.bam")
                output_directory = os.path.join("metrics", sample.name)
                
                job = concat_jobs([
                        Job(command="mkdir -p " + output_directory, removable_files=[output_directory]),
                        picard.collect_multiple_metrics(alignment_file, os.path.join(output_directory,sample.name),reference_file),
                        picard.collect_rna_metrics(alignment_file, os.path.join(output_directory,sample.name+".picard_rna_metrics"))
                ],name="picard_rna_metrics."+ sample.name)
                jobs.append(job)
        
        return jobs

    def bam_hard_clip(self):
        """
        Generate a hardclipped version of the bam for the toxedo suite which doesn't support this official sam feature.
        """
        
        jobs = []
        for sample in self.samples:
            alignment_input = os.path.join("alignment", sample.name, sample.name + ".sorted.bam")
            alignment_output = os.path.join("alignment", sample.name, sample.name + ".sorted.hardClip.bam")
            job=pipe_jobs([
                samtools.view(
                    alignment_input,
                    None,
                    "-h"
                ),
		Job(
                    [None],
                    [alignment_output],
                    # awk to transform soft clip into hard clip for tuxedo suite
                    command="""\
awk 'BEGIN {{OFS="\\t"}} {{if (substr($1,1,1)=="@") {{print;next}}; split($6,C,/[0-9]*/); split($6,L,/[SMDIN]/); if (C[2]=="S") {{$10=substr($10,L[1]+1); $11=substr($11,L[1]+1)}}; if (C[length(C)]=="S") {{L1=length($10)-L[length(L)-1]; $10=substr($10,1,L1); $11=substr($11,1,L1); }}; gsub(/[0-9]*S/,"",$6); print}}' """.format()
                ),
                samtools.view(
                    "-",
                    alignment_output,
                    "-hbS"
                ),
            ])
            job.name="tuxedo_hard_clip."+ sample.name
            jobs.append(job)
        return jobs    

    def cufflinks(self):
        """        Compute RNA-Seq data expression using [cufflinks](http://cole-trapnell-lab.github.io/cufflinks/cufflinks/).
        Warning: It needs to use a hard clipped bam file while Tuxedo tools do not support official soft clip SAM format
        """

        jobs = []
        
        gtf = config.param('cufflinks','gtf', type='filepath')
        for sample in self.samples:
            input_bam = os.path.join("alignment", sample.name, sample.name + ".sorted.hardClip.bam")
            output_directory = os.path.join("cufflinks", sample.name)

            # De Novo FPKM
            job = cufflinks.cufflinks(input_bam, output_directory, gtf)
            job.removable_files = ["cufflinks"]
            job.name = "cufflinks."+sample.name
            jobs.append(job)

        return jobs

    def cuffmerge(self):
        """
        Merge assemblies into a master transcriptome reference using [cuffmerge](http://cole-trapnell-lab.github.io/cufflinks/cuffmerge/).
        """
        output_directory = os.path.join("cufflinks", "AllSamples")
        sample_file = os.path.join("cufflinks", "cuffmerge.samples.txt")
        input_gtfs = [os.path.join("cufflinks", sample.name, "transcripts.gtf") for sample in self.samples]
        gtf = config.param('cuffmerge','gtf', type='filepath')
        
        
        job = concat_jobs([
            Job(command="mkdir -p " + output_directory),
            Job(input_gtfs, [sample_file], command="""\
`cat > {sample_file} << END
{sample_rows}
END
  
`""".format(sample_rows="\n".join(input_gtfs), sample_file=sample_file)),
            cufflinks.cuffmerge(sample_file, output_directory, gtf_file=gtf)],
            name="cuffmerge")
        
        return [job]

    def cuffquant(self):
        """
        Compute expression profiles (abundances.cxb) using [cuffquant](http://cole-trapnell-lab.github.io/cufflinks/cuffquant/).
        Warning: It needs to use a hard clipped bam file while Tuxedo tools do not support official soft clip SAM format
        """

        jobs = []
        
        gtf = os.path.join("cufflinks", "AllSamples","merged.gtf")
        
        for sample in self.samples:
            input_bam = os.path.join("alignment", sample.name, sample.name + ".sorted.hardClip.bam")
            output_directory = os.path.join("cufflinks", sample.name)

            #Quantification
            job = cufflinks.cuffquant(input_bam, output_directory, gtf)
            job.name = "cuffquant."+sample.name
            jobs.append(job)

        return jobs

    def cuffdiff(self):
        """
        [Cuffdiff](http://cole-trapnell-lab.github.io/cufflinks/cuffdiff/) is used to calculate differential transcript expression levels and test them for significant differences.
        """

        jobs = []

        fpkm_directory = "cufflinks"
        gtf = os.path.join(fpkm_directory, "AllSamples","merged.gtf")


        # Perform cuffdiff on each design contrast
        for contrast in self.contrasts:
            job = cufflinks.cuffdiff(
                # Cuffdiff input is a list of lists of replicate bams per control and per treatment
                [[os.path.join(fpkm_directory, sample.name, "abundances.cxb") for sample in group] for group in contrast.controls, contrast.treatments],
                gtf,
                os.path.join("cuffdiff", contrast.name)
            )
            job.removable_files = ["cuffdiff"]
            job.name = "cuffdiff." + contrast.name
            jobs.append(job)

        return jobs

    def cuffnorm(self):
        """        Global normalization of RNA-Seq expression levels using [Cuffnorm](http://cole-trapnell-lab.github.io/cufflinks/cuffnorm/).
        """

        jobs = []

        fpkm_directory = "cufflinks"
        gtf = os.path.join(fpkm_directory, "AllSamples","merged.gtf")
        sample_labels = ",".join([sample.name for sample in self.samples])

        # Perform cuffnorm using every samples
        job = cufflinks.cuffnorm([os.path.join(fpkm_directory, sample.name, "abundances.cxb") for sample in self.samples],
             gtf,
             "cuffnorm",sample_labels)
        job.removable_files = ["cuffnorm"]
        job.name = "cuffnorm" 
        jobs.append(job)

        report_file = os.path.join("report", "RNAseq_AS.cuffdiff.md")
        jobs.append(
            Job(
                [os.path.join(fpkm_directory, sample.name, "abundances.cxb") for sample in self.samples],
                [report_file],
                [['bedtools', 'module_pandoc']],
                command="""\
mkdir -p report && \\
cp cuffdiff/Contrast1/splicing.diff report && \\
pandoc --to=markdown \\
--template {report_template_dir}/{basename_report_file} \\
{report_template_dir}/{basename_report_file} \\
> {report_file}""".format(
                    report_template_dir=self.report_template_dir,
                    basename_report_file=os.path.basename(report_file),
                    report_file=report_file
                    ),
                report_files=[report_file],
                name="cuffdiff_report")
            )

        return jobs

    def miso_index(self):
        """
        Creates the indexed folder needed to compute the Psi values from a gff file

        """
        jobs = []
        
        gffs = ['A3SS.hg19.gff3','A5SS.hg19.gff3','MXE.hg19.gff3','RI.hg19.gff3','SE.hg19.gff3']

        for gff in gffs:
            
            indexed_folder = 'miso/indexed_' + gff.split('.')[0]
            input_gff_file = os.path.join('/hpf/projects/brudno/vincent/data/hg19', gff)

            job = miso_funcs.make_index(indexed_folder, input_gff_file)
            job.name = "miso_index." + gff.split('.')[0]

            jobs.append(job)

        return jobs

    def miso_paired_end(self):
        """
        Calculates the mean insert length and standard deviation for paired end reads
        """

        jobs = []

        if config.param('DEFAULT', 'library_type', type='string', required=True) == 'paired':
            for sample in self.samples:

                alignment_file = os.path.join('alignment', sample.name, sample.name + '.sorted.bam')

                job = miso_funcs.paired_end(alignment_file, sample.name)

                job.name = 'miso_paired_end.' + sample.name
                jobs.append(job)

        return jobs

    def miso_psi(self):
        """
        FIrst it uses SAMtools to create the header file for the sorted bams. Then it computes the Psi values of the bams.

        """

        jobs = []

        for sample in self.samples:

            gffs = ['A3SS','A5SS','MXE','RI','SE']
            for gff in gffs:

                indexed_folder = 'miso/indexed_' + gff
                alignment_file_name = os.path.join("alignment", sample.name, sample.name + ".sorted.bam")
                bam_header = alignment_file_name + ".bai"
                output_folder_name = os.path.join("miso", sample.name, gff)

                job = miso_funcs.compute_psi(alignment_file_name, indexed_folder, bam_header, output_folder_name, sample.name)
                job.name = 'miso_psi.' + sample.name + '.' + gff
                jobs.append(job)

        return jobs

    def miso_summarize(self):
       
        jobs = []

        for sample in self.samples:

            gffs = ['A3SS','A5SS','MXE','RI','SE']
            for gff in gffs:

                summary_folder_prefix = os.path.join("miso", sample.name, gff)

                job = miso_funcs.summary(summary_folder_prefix, os.path.join(summary_folder_prefix, "summary", sample.name + ".miso_summary"))
                job.name = "miso_summarize." + sample.name + '.' + gff
                jobs.append(job)

        return jobs

    def miso_diff(self):
        jobs = []

        summaries = []
        for sample in self.samples:

            gffs = ['A3SS','A5SS','MXE','RI','SE']
            for gff in gffs:

                summaries.append(os.path.join('miso', sample.name, gff, 'summary', sample.name + '.miso_summary'))

        for contrast in self.contrasts:

            control_sample = (contrast.controls[0]).name
            treatment_sample = (contrast.treatments[0]).name
            output_folder = os.path.join("miso", "comparisons", contrast.name)
            output = os.path.join(output_folder, control_sample + "_vs_" + treatment_sample, "bayes-factors", control_sample + "_vs_" + treatment_sample + ".miso_bf")

            job = miso_funcs.compare(
                control_sample,
                treatment_sample,
                output_folder,
                output,
                summaries)

            job.name = "miso_diff." + contrast.name
            jobs.append(job)

        return jobs

    def miso_plot(self):

        jobs = []

        report_dependencies = []
        event_plots = ''

        plot_type = config.param('miso_plot', 'plot_type', type='string', required=True)

        sample_names = '['
        bam_files = '['
        for sample in self.samples:
            sample_names += '\"' +  sample.name + '\", ' 
            bam_files += '\"alignment/' + sample.name + '/' + sample.name + '.sorted.bam\", '
        sample_names = sample_names[:-2] + ']'
        bam_files = bam_files[:-2] + ']'
        job = miso_funcs.plot_settings(sample_names, bam_files)
        job.name = 'miso_plot_settings'
        jobs.append(job)

        if '--plot-insert-len' in plot_type:

            if config.param('DEFAULT', 'library_type', type='string', required=True) == 'paired':

                for sample in self.samples:
                    insert_file = os.path.join('miso', sample.name, 'insert-dist', sample.name + '.sorted.bam.insert_len')

                    job = miso_funcs.plot_insert_len(insert_file)
                    job.name = 'miso_plot.insert_len.' + sample.name
                    jobs.append(job)

        if '--plot-bf-dist' in plot_type:

            for contrast in self.contrasts:

                    control_sample = (contrast.controls[0]).name
                    treatment_sample = (contrast.treatments[0]).name
                    bf_file = os.path.join('miso', 'comparisons', contrast.name, control_sample + '_vs_' + treatment_sample, 'bayes-factors', control_sample + '_vs_' + '.miso_bf')

                    job = miso_funcs.plot_bf_dist(bf_file)
                    job.name = 'miso_plot.bf_dist.' + contrast.name
                    jobs.append(job)

        if '--plot-event' in plot_type:

            events_names = config.param('miso_plot', 'events_names', type='string', required=True).split()

            dependencies = []
            for sample in self.samples:
                
                gffs = ['A3SS','A5SS','MXE','RI','SE']
                for gff in gffs:
                    dependencies.append( os.path.join('miso', sample.name, gff, 'summary', sample.name + '.miso_summary'))
            
            i = 0
            
            for event_name in events_names:
            
                i += 1
                report_dependencies.append(os.path.join('miso', 'plots', event_name + '.pdf'))
                job = miso_funcs.plot_event(dependencies, event_name)
                job.name = 'miso_plot.plot_event.' + str(i)
                jobs.append(job)

################################# report generation

        report_file = os.path.join("report", "RNAseq_AS.miso_plot.md")
        report_template_dir = self.report_template_dir
        basename_report_file = os.path.basename(report_file)

        for contrast in self.contrasts:
            control_sample = (contrast.controls[0]).name
            treatment_sample = (contrast.treatments[0]).name
            output_bf = os.path.join("miso", "comparisons", contrast.name, control_sample + "_vs_" + treatment_sample, "bayes-factors", control_sample + "_vs_" + treatment_sample + ".miso_bf")
            report_dependencies.append(output_bf)

        job = miso_funcs.report(
            report_dependencies,
            report_file,
            report_template_dir,
            basename_report_file
            )

        jobs.append(job)

        return jobs

    def bam_to_fastq(self):
        jobs = []

        for sample in self.samples:
            alignment_file_prefix = os.path.join("alignment", sample.name, sample.name)

            job = vtools.make_fastq(
                alignment_file_prefix + ".sorted.bam",
                alignment_file_prefix + ".fastq"
                )

            job.name = "bam_to_fastq." + sample.name
            jobs.append(job)

        return jobs

    def vast_tools_align(self):
        jobs = []

        for sample in self.samples:
            input_fastq = os.path.join("alignment", sample.name, sample.name + ".fastq")
            align_out = os.path.join("vast_out", "to_combine", sample.name + ".eej2")

            job = vtools.alignment(input_fastq, align_out)
            
            job.name = "vast_tools_align." + sample.name
            jobs.append(job)

        return jobs

    def vast_tools_combine(self):
        jobs = []
        
        align_out_list = []
        for sample in self.samples:
            align_out_list.append(os.path.join("vast_out", "to_combine", sample.name + ".eej2"))

        species = config.param('vast_tools_align', 'species', type='string', required=True)
        full_inclusion = "vast_out/INCLUSION_LEVELS_FULL-" + species + str(len(align_out_list)) + ".tab"
        
        job = vtools.comb(align_out_list, full_inclusion)

        job.name = "vast_tools_combine"
        jobs.append(job)

        return jobs

    def vast_tools_diff(self):
        jobs = []
                
        num_samples = 0
        for sample in self.samples:
                num_samples += 1

        species = config.param('vast_tools_align', 'species', type='string', required=True)
        full_inclusion = "vast_out/INCLUSION_LEVELS_FULL-" + species + str(num_samples) + ".tab"
       
        for contrast in self.contrasts:

            control_samples_names = []
            for control in contrast.controls:
                control_samples_names.append(control.name)

            treatment_samples_names = []
            for treatment in contrast.treatments:
                treatment_samples_names.append(treatment.name)

            filtered_inclusion_table = os.path.join("vast_out", "INCLUSION-FILTERED-" + contrast.name + ".tab" )

            job = vtools.differential_splicing(control_samples_names, treatment_samples_names, full_inclusion, filtered_inclusion_table)

            job.name = "vast_tools_diff." + contrast.name
            jobs.append(job)
        
        return jobs

    def vast_tools_plot(self):
        
        jobs = []
        
        num_samples = 0
        for sample in self.samples:
                num_samples += 1

        species = config.param('vast_tools_align', 'species', type='string', required=True)
        full_inclusion = 'vast_out/INCLUSION_LEVELS_FULL-' + species + str(num_samples) + '.tab'
        events_list = config.param('vast_tools_plot', 'significant_events', type='string', required=True).split()

        job = vtools.plots(full_inclusion, events_list)
        job.name = 'vast_tools_plot'
        jobs.append(job)

################### report gen

        report_file = os.path.join('report', 'RNAseq_AS.vtools_plot.md')
        report_template_dir = self.report_template_dir
        basename_report_file = os.path.basename(report_file)

        job = vtools.report(
            report_file,
            report_template_dir,
            basename_report_file,
            full_inclusion,
            )
        
        jobs.append(job)

        return jobs

    def jctseq_raw_counts(self):

        jobs = []

        for sample in self.samples:
            alignment_file = os.path.join("alignment", sample.name, sample.name + ".sorted.bam")
            raw_counts_dir = os.path.join("jctseq", "rawCts", sample.name)
            job = jctseq.raw_counts(
                alignment_file,
                raw_counts_dir,
                )

            job.name = "jctseq_raw_counts." + sample.name
            jobs.append(job)

        return jobs

    def jctseq_make_gff(self):

        jobs = []

        job = jctseq.make_gff("jctseq/jctseq.gff.gz")
        job.name = "jctseq_make_gff"
        jobs.append(job)

        return jobs

    def jctseq_diff_prep(self):
        
        jobs = []
        gff = os.path.join('jctseq', 'jctseq.gff.gz')
        
        for contrast in self.contrasts:
            folder = os.path.join('jctseq', contrast.name)

            job = jctseq.diff_prep(folder, contrast.name, gff)
            job.name = 'jctseq_diff_prep.' + contrast.name
            jobs.append(job)

        return jobs

    def jctseq_diff(self):

        jobs = []

        raw_counts = []
        for sample in self.samples:
            raw_counts.append(os.path.join('jctseq', 'rawCts', sample.name, 'QC.QORTS_COMPLETED_OK'))

        report_dependencies = []
        for contrast in self.contrasts:
            
            jscs_file = os.path.join('jctseq', contrast.name, 'jscs1.r')
            report_dependencies.append(os.path.join('jctseq', contrast.name, 'plots'))

            job = jctseq.jscs_diff(
                'jctseq/jctseq.gff.gz', 
                os.path.join('jctseq', contrast.name, 'plots'),
                os.path.join('jctseq', contrast.name, 'results'),
                jscs_file,
                raw_counts)

            job.name = 'jctseq_diff.' + contrast.name
            jobs.append(job)

###### report
        report_file = os.path.join('report', 'RNAseq_AS.jctseq_diff.md')
        report_template_dir = self.report_template_dir
        basename_report_file = os.path.basename(report_file)

        job = jctseq.report(
            report_dependencies,
            report_file,
            report_template_dir,
            basename_report_file
            )
        
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
	    self.bam_hard_clip,
            self.cufflinks,
            self.cuffmerge,
            self.cuffquant,
            self.cuffdiff,
            self.cuffnorm,
	    self.miso_index,
            self.miso_paired_end,
            self.miso_psi,
            self.miso_summarize,
            self.miso_diff,
            self.miso_plot,
	    self.bam_to_fastq,
            self.vast_tools_align,
            self.vast_tools_combine,
            self.vast_tools_diff,
            self.vast_tools_plot,
            self.jctseq_raw_counts,
            self.jctseq_make_gff,
            self.jctseq_diff_prep,
            self.jctseq_diff
            ]

if __name__ == '__main__':
    RnaSeq()
