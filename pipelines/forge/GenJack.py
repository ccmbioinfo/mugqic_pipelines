#!/usr/bin/env python

# Python Standard Modules
import logging
import math
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
from bfx.sample_info import *
from bfx.sequence_dictionary import *

# These tools will change for this pipeline
from bfx import trimmomatic
from bfx import fast
from bfx import bwa
from bfx import gatk
from bfx import gq_seq_utils
from bfx import igvtools
from bfx import metrics
from bfx import picard
from bfx import samtools
from bfx import tools
from bfx import vcftools
from bfx import forge_tools
from bfx import clean
from bfx import cramtools
from bfx import bvatools
from bfx import bcftools
from bfx import tabix
from bfx import annovar
from bfx import vt
from bfx import vep
from bfx import gemini

from pipelines import common

log = logging.getLogger(__name__)

class Forge(common.Illumina):
    """
    Forge Exome Sequencing Pipeline
    ===============================

    """
    def __init__(self):
        # Add pipeline specific arguments
        self._lastPGStep = {}
        super(Forge, self).__init__()


    @property
    def sample_info(self):
        if not hasattr(self, "_sample_info"):
            self._sample_info = {}

            # Create a dictionary of dictionaries holding all the option values
            sample_info_file = config.param("DEFAULT", "sample_info_file", type="filepath")
            with open(sample_info_file, "r") as f:
                for line in f.readlines():
                    line = line.strip()
                    # Ignore comments
                    if not line.startswith("#"):
                        split_sample = line.split("\t")
                        # Make sure the format is correct
                        if len(split_sample) > 1:
                            param = split_sample[1].split("=")
                            if split_sample[0] in self._sample_info.keys():
                                self._sample_info[split_sample[0]][param[0]] = param[1] if len(param) > 1 else 1 # If not equal to anything, set to true b/c calle
                            else:
                                self._sample_info[split_sample[0]] = {param[0]: param[1]} if len(param) > 1 else {param[0]: 1}

        return self._sample_info


    @property
    def sequence_dictionary(self):
        if not hasattr(self, "_sequence_dictionary"):
            self._sequence_dictionary = parse_sequence_dictionary_file(config.param("DEFAULT", "genome_dictionary", type="filepath"))
        return self._sequence_dictionary


    def get_ver(self, module):
        return config.param("DEFAULT", "module_" + module).split("/")[-1]


    def trimmomatic(self):
        """
        Copied from common.py - need to add modification of saving the command to the PG file

        Raw reads quality trimming and removing of Illumina adapters is performed using [Trimmomatic](http://www.usadellab.org/cms/index.php?page=trimmomatic).
        If an adapter FASTA file is specified in the config file (section 'trimmomatic', param 'adapter_fasta'),
        it is used first. Else, 'Adapter1' and 'Adapter2' columns from the readset file are used to create
        an adapter FASTA file, given then to Trimmomatic. For PAIRED_END readsets, readset adapters are
        reversed-complemented and swapped, to match Trimmomatic Palindrome strategy. For SINGLE_END readsets,
        only Adapter1 is used and left unchanged.

        This step takes as input files:

        1. FASTQ files from the readset file if available
        2. Else, FASTQ output files from previous picard_sam_to_fastq conversion of BAM files
        """
        jobs = []
        self._lastPGStep["default"] = "ref_genome_indexing" # For PG File

        for readset in self.readsets:
            trim_directory = os.path.join("trim", readset.sample.name)
            trim_file_prefix = os.path.join(trim_directory, readset.name + ".trim.")
            trim_log = trim_file_prefix + "log"

            # Use adapter FASTA in config file if any, else create it from readset file
            adapter_fasta = config.param("trimmomatic", "adapter_fasta", required=False, type="filepath")
            adapter_job = None
            if not adapter_fasta:
                adapter_fasta = trim_file_prefix + "adapters.fa"
                if readset.run_type == "PAIRED_END":
                    if readset.adapter1 and readset.adapter2:
                        # WARNING: Reverse-complement and swap readset adapters for Trimmomatic Palindrome strategy
                        adapter_job = Job(command="""\
`cat > {adapter_fasta} << END
>Prefix/1
{sequence1}
>Prefix/2
{sequence2}
END
`""".format(adapter_fasta=adapter_fasta, sequence1=readset.adapter2.translate(string.maketrans("ACGTacgt","TGCAtgca"))[::-1], sequence2=readset.adapter1.translate(string.maketrans("ACGTacgt","TGCAtgca"))[::-1]))
                    else:
                        raise Exception("Error: missing adapter1 and/or adapter2 for PAIRED_END readset \"" + readset.name + "\", or missing adapter_fasta parameter in config file!")
                elif readset.run_type == "SINGLE_END":
                    if readset.adapter1:
                        adapter_job = Job(command="""\
`cat > {adapter_fasta} << END
>Single
{sequence}
END
`""".format(adapter_fasta=adapter_fasta, sequence=readset.adapter1))
                    else:
                        raise Exception("Error: missing adapter1 for SINGLE_END readset \"" + readset.name + "\", or missing adapter_fasta parameter in config file!")


            trim_stats = trim_file_prefix + "stats.csv"
            if readset.run_type == "PAIRED_END":
                candidate_input_files = [[readset.fastq1, readset.fastq2]]
                if readset.bam:
                    candidate_input_files.append([re.sub("\.bam$", ".pair1.fastq.gz", readset.bam), re.sub("\.bam$", ".pair2.fastq.gz", readset.bam)])
                [fastq1, fastq2] = self.select_input_files(candidate_input_files)
                job = trimmomatic.trimmomatic(
                    fastq1,
                    fastq2,
                    trim_file_prefix + "pair1.fastq.gz",
                    trim_file_prefix + "single1.fastq.gz",
                    trim_file_prefix + "pair2.fastq.gz",
                    trim_file_prefix + "single2.fastq.gz",
                    None,
                    readset.quality_offset,
                    adapter_fasta,
                    trim_log
                )
            elif readset.run_type == "SINGLE_END":
                candidate_input_files = [[readset.fastq1]]
                if readset.bam:
                    candidate_input_files.append([re.sub("\.bam$", ".single.fastq.gz", readset.bam)])
                [fastq1] = self.select_input_files(candidate_input_files)
                job = trimmomatic.trimmomatic(
                    fastq1,
                    None,
                    None,
                    None,
                    None,
                    None,
                    trim_file_prefix + "single.fastq.gz",
                    readset.quality_offset,
                    adapter_fasta,
                    trim_log
                )
            else:
                raise Exception("Error: run type \"" + readset.run_type +
                "\" is invalid for readset \"" + readset.name + "\" (should be PAIRED_END or SINGLE_END)!")

            trim_job = job # Keep the original trimmomatic job for the command to be passed to the PG file
            trim_ver = config.param("DEFAULT", "module_trimmomatic").split("/")[-1]

            if adapter_job:
                job = concat_jobs([adapter_job, job])

            jobs.append(concat_jobs([
                # Trimmomatic does not create output directory by default
                Job(command="mkdir -p " + trim_directory),
                job,
                forge_tools.add_pg_file(trim_directory + "/" + readset.sample.name + ".pg.txt", "trim-" + readset.name, trim_ver, "\"" + trim_job.command + "\"", self._lastPGStep["default"])
            ], name="trimmomatic." + readset.name))

            self._lastPGStep[readset.name] = "trim-" + readset.name

        return jobs


##############################################################################################################################

    @property
    def sequence_dictionary(self):
        if not hasattr(self, "_sequence_dictionary"):
            self._sequence_dictionary = parse_sequence_dictionary_file(config.param('DEFAULT', 'genome_dictionary', type='filepath'))
        return self._sequence_dictionary


##############################################################################################################################


    def bwa_mem_picard_sort_sam(self):
        """
        The filtered reads are aligned to a reference genome. The alignment is done per sequencing readset.
        The alignment software used is [BWA](http://bio-bwa.sourceforge.net/) with algorithm: bwa mem.
        BWA output BAM files are then sorted by coordinate using [Picard](http://broadinstitute.github.io/picard/).

        This step takes as input files:

        1. Trimmed FASTQ files if available
        2. Else, FASTQ files from the readset file if available
        3. Else, FASTQ output files from previous picard_sam_to_fastq conversion of BAM files
        """

        jobs = []
        for readset in self.readsets:
            trim_file_prefix = os.path.join("trim", readset.sample.name, readset.name + ".trim.")
            alignment_directory = os.path.join("alignment", readset.sample.name)
            readset_bam = os.path.join(alignment_directory, readset.name, readset.name + ".sorted.bam")

            # Find input readset FASTQs first from previous trimmomatic job, then from original FASTQs in the readset sheet
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



            job = concat_jobs([
                Job(command="mkdir -p " + os.path.dirname(readset_bam)),
                pipe_jobs([
                    bwa.mem(
                        fastq1,
                        fastq2,
                        read_group="'@RG" + \
                            "\tID:" + readset.name + \
                            "\tSM:" + readset.sample.name + \
                            ("\tLB:" + readset.library if readset.library else "") + \
                            ("\tPU:run" + readset.run + "_" + readset.lane if readset.run and readset.lane else "") + \
                            ("\tCN:" + config.param('bwa_mem', 'sequencing_center') if config.param('bwa_mem', 'sequencing_center', required=False) else "") + \
                            "\tPL:Illumina" + \
                            "'"
                    ),
                    picard.sort_sam(
                        "/dev/stdin",
                        readset_bam,
                        "coordinate"
                    )
                ])
            ], name="bwa_mem_picard_sort_sam." + readset.name)

            jobs.append(job)

        report_file = os.path.join("report", "DnaSeq.bwa_mem_picard_sort_sam.md")
        jobs.append(
            Job(
                [os.path.join("alignment", readset.sample.name, readset.name, readset.name + ".sorted.bam") for readset in self.readsets],
                [report_file],
                [['bwa_mem_picard_sort_sam', 'module_pandoc']],
                command="""\
mkdir -p report && \\
pandoc --to=markdown \\
  --template {report_template_dir}/{basename_report_file} \\
  --variable scientific_name="{scientific_name}" \\
  --variable assembly="{assembly}" \\
  {report_template_dir}/{basename_report_file} \\
  > {report_file}""".format(
                    scientific_name=config.param('bwa_mem_picard_sort_sam', 'scientific_name'),
                    assembly=config.param('bwa_mem_picard_sort_sam', 'assembly'),
                    report_template_dir=self.report_template_dir,
                    basename_report_file=os.path.basename(report_file),
                    report_file=report_file
                ),
                report_files=[report_file],
                name="bwa_mem_picard_sort_sam_report")
        )

        return jobs



##############################################################################################################################



    def picard_merge_sam_files(self):
        """
        BAM readset files are merged into one file per sample. Merge is done using [Picard](http://broadinstitute.github.io/picard/).

        This step takes as input files:

        1. Aligned and sorted BAM output files from previous bwa_mem_picard_sort_sam step if available
        2. Else, BAM files from the readset file
        """

        jobs = []
        for sample in self.samples:
            alignment_directory = os.path.join("alignment", sample.name)
            # Find input readset BAMs first from previous bwa_mem_picard_sort_sam job, then from original BAMs in the readset sheet.
            readset_bams = self.select_input_files([[os.path.join(alignment_directory, readset.name, readset.name + ".sorted.bam") for readset in sample.readsets], [readset.bam for readset in sample.readsets]])
            sample_bam = os.path.join(alignment_directory, sample.name + ".sorted.bam")

            mkdir_job = Job(command="mkdir -p " + os.path.dirname(sample_bam))

            # If this sample has one readset only, create a sample BAM symlink to the readset BAM, along with its index.
            if len(sample.readsets) == 1:
                readset_bam = readset_bams[0]
                if os.path.isabs(readset_bam):
                    target_readset_bam = readset_bam
                else:
                    target_readset_bam = os.path.relpath(readset_bam, alignment_directory)
                readset_index = re.sub("\.bam$", ".bai", readset_bam)
                target_readset_index = re.sub("\.bam$", ".bai", target_readset_bam)
                sample_index = re.sub("\.bam$", ".bai", sample_bam)

                job = concat_jobs([
                    mkdir_job,
                    Job([readset_bam], [sample_bam], command="ln -s -f " + target_readset_bam + " " + sample_bam, removable_files=[sample_bam]),
                    Job([readset_index], [sample_index], command="ln -s -f " + target_readset_index + " " + sample_index, removable_files=[sample_index])
                ], name="symlink_readset_sample_bam." + sample.name)

            elif len(sample.readsets) > 1:
                job = concat_jobs([
                    mkdir_job,
                    picard.merge_sam_files(readset_bams, sample_bam)
                ])
                job.name = "picard_merge_sam_files." + sample.name

            jobs.append(job)

        return jobs

##############################################################################################################################


    def gatk_indel_realigner(self):
        """
        Insertion and deletion realignment is performed on regions where multiple base mismatches
        are preferred over indels by the aligner since it can appear to be less costly by the algorithm.
        Such regions will introduce false positive variant calls which may be filtered out by realigning
        those regions properly. Realignment is done using [GATK](https://www.broadinstitute.org/gatk/).
        The reference genome is divided by a number regions given by the `nb_jobs` parameter.
        """

        jobs = []

        nb_jobs = config.param('gatk_indel_realigner', 'nb_jobs', type='posint')
        if nb_jobs > 50:
            log.warning("Number of realign jobs is > 50. This is usually much. Anything beyond 20 can be problematic.")

        for sample in self.samples:
            alignment_directory = os.path.join("alignment", sample.name)
            realign_directory = os.path.join(alignment_directory, "realign")
            input = os.path.join(alignment_directory, sample.name + ".sorted.bam")

            if nb_jobs == 1:
                realign_prefix = os.path.join(realign_directory, "all")
                realign_intervals = realign_prefix + ".intervals"
                output_bam = realign_prefix + ".bam"
                sample_output_bam = os.path.join(alignment_directory, sample.name + ".realigned.qsorted.bam")
                jobs.append(concat_jobs([
                    Job(command="mkdir -p " + realign_directory, removable_files=[realign_directory]),
                    gatk.realigner_target_creator(input, realign_intervals),
                    gatk.indel_realigner(input, output_bam, target_intervals=realign_intervals),
                    # Create sample realign symlink since no merging is required
                    Job([output_bam], [sample_output_bam], command="ln -s -f " + os.path.relpath(output_bam, os.path.dirname(sample_output_bam)) + " " + sample_output_bam)
                ], name="gatk_indel_realigner." + sample.name))

            else:
                # The first sequences are the longest to process.
                # Each of them must be processed in a separate job.
                unique_sequences_per_job = [sequence['name'] for sequence in self.sequence_dictionary[0:min(nb_jobs - 1, len(self.sequence_dictionary))]]

                # Create one separate job for each of the first sequences
                for sequence in unique_sequences_per_job:
                    realign_prefix = os.path.join(realign_directory, sequence)
                    realign_intervals = realign_prefix + ".intervals"
                    intervals=[sequence]
                    if unique_sequences_per_job.index(sequence) == 0:
                        intervals.append("unmapped")
                    output_bam = realign_prefix + ".bam"
                    jobs.append(concat_jobs([
                        # Create output directory since it is not done by default by GATK tools
                        Job(command="mkdir -p " + realign_directory, removable_files=[realign_directory]),
                        gatk.realigner_target_creator(input, realign_intervals, intervals=[sequence]),
                        gatk.indel_realigner(input, output_bam, target_intervals=realign_intervals, intervals=intervals)
                    ], name="gatk_indel_realigner." + sample.name + "." + sequence))

                # Create one last job to process the last remaining sequences and 'others' sequences
                realign_prefix = os.path.join(realign_directory, "others")
                realign_intervals = realign_prefix + ".intervals"
                output_bam = realign_prefix + ".bam"
                jobs.append(concat_jobs([
                    # Create output directory since it is not done by default by GATK tools
                    Job(command="mkdir -p " + realign_directory, removable_files=[realign_directory]),
                    gatk.realigner_target_creator(input, realign_intervals, exclude_intervals=unique_sequences_per_job),
                    gatk.indel_realigner(input, output_bam, target_intervals=realign_intervals, exclude_intervals=unique_sequences_per_job)
                ], name="gatk_indel_realigner." + sample.name + ".others"))

        return jobs


#####################################################################################################################
    def merge_realigned(self):
        """
        BAM files of regions of realigned reads are merged per sample using [Picard](http://broadinstitute.github.io/picard/).
        """

        jobs = []

        nb_jobs = config.param('gatk_indel_realigner', 'nb_jobs', type='posint')

        for sample in self.samples:
            alignment_directory = os.path.join("alignment", sample.name)
            realign_directory = os.path.join(alignment_directory, "realign")
            merged_realigned_bam = os.path.join(alignment_directory, sample.name + ".realigned.qsorted.bam")

            # if nb_jobs == 1, symlink has been created in indel_realigner and merging is not necessary
            if nb_jobs > 1:
                realigned_bams = [os.path.join(realign_directory, sequence['name'] + ".bam") for sequence in self.sequence_dictionary[0:min(nb_jobs - 1, len(self.sequence_dictionary))]]
                realigned_bams.append(os.path.join(realign_directory, "others.bam"))

                job = picard.merge_sam_files(realigned_bams, merged_realigned_bam)
                job.name = "merge_realigned." + sample.name
                jobs.append(job)

        report_file = os.path.join("report", "DnaSeq.gatk_indel_realigner.md")
        jobs.append(
            Job(
                [os.path.join("alignment", sample.name, sample.name + ".realigned.qsorted.bam") for sample in self.samples],
                [report_file],
                command="""\
mkdir -p report && \\
cp \\
  {report_template_dir}/{basename_report_file} \\
  {report_file}""".format(
                    report_template_dir=self.report_template_dir,
                    basename_report_file=os.path.basename(report_file),
                    report_file=report_file
                ),
                report_files=[report_file],
                name="merge_realigned_report")
        )

        return jobs

############################################################################################################################

    def fix_mate_by_coordinate(self):
        """
        Fix the read mates. Once local regions are realigned, the read mate coordinates of the aligned reads
        need to be recalculated since the reads are realigned at positions that differ from their original alignment.
        Fixing the read mate positions is done using [BVATools](https://bitbucket.org/mugqic/bvatools).
        """

        jobs = []
        for sample in self.samples:
            alignment_file_prefix = os.path.join("alignment", sample.name, sample.name + ".")
            input = alignment_file_prefix + "realigned.qsorted.bam"
            output_prefix = alignment_file_prefix + "matefixed.sorted"
            jobs.append(concat_jobs([
                bvatools.groupfixmate(input, output_prefix + ".tmp.bam"),
                picard.sort_sam(output_prefix + ".tmp.bam", output_prefix+ ".bam"),
            ], name="fix_mate_by_coordinate." + sample.name))

        report_file = os.path.join("report", "DnaSeq.fix_mate_by_coordinate.md")
        jobs.append(
            Job(
                [os.path.join("alignment", sample.name, sample.name + ".matefixed.sorted.bam") for sample in self.samples],
                [report_file],
                command="""\
mkdir -p report && \\
cp \\
  {report_template_dir}/{basename_report_file} \\
  {report_file}""".format(
                    report_template_dir=self.report_template_dir,
                    basename_report_file=os.path.basename(report_file),
                    report_file=report_file
                ),
                report_files=[report_file],
                name="fix_mate_by_coordinate_report")
        )

        return jobs



############################################################################################################################

    def picard_mark_duplicates(self):
        """
        Mark duplicates. Aligned reads per sample are duplicates if they have the same 5' alignment positions
        (for both mates in the case of paired-end reads). All but the best pair (based on alignment score)
        will be marked as a duplicate in the BAM file. Marking duplicates is done using [Picard](http://broadinstitute.github.io/picard/).
        """

        jobs = []
        for sample in self.samples:
            alignment_file_prefix = os.path.join("alignment", sample.name, sample.name + ".")
            input = alignment_file_prefix + "matefixed.sorted.bam"
            output = alignment_file_prefix + "sorted.dup.bam"
            metrics_file = alignment_file_prefix + "sorted.dup.metrics"

            job = picard.mark_duplicates([input], output, metrics_file)
            job.name = "picard_mark_duplicates." + sample.name
            jobs.append(job)

        report_file = os.path.join("report", "DnaSeq.picard_mark_duplicates.md")
        jobs.append(
            Job(
                [os.path.join("alignment", sample.name, sample.name + ".sorted.dup.bam") for sample in self.samples],
                [report_file],
                command="""\
mkdir -p report && \\
cp \\
  {report_template_dir}/{basename_report_file} \\
  {report_file}""".format(
                    report_template_dir=self.report_template_dir,
                    basename_report_file=os.path.basename(report_file),
                    report_file=report_file
                ),
                report_files=[report_file],
                name="picard_mark_duplicates_report")
        )

        return jobs

##########################################################################################################################

    def recalibration(self):
        """
        Recalibrate base quality scores of sequencing-by-synthesis reads in an aligned BAM file. After recalibration,
        the quality scores in the QUAL field in each read in the output BAM are more accurate in that
        the reported quality score is closer to its actual probability of mismatching the reference genome.
        Moreover, the recalibration tool attempts to correct for variation in quality with machine cycle
        and sequence context, and by doing so, provides not only more accurate quality scores but also
        more widely dispersed ones.
        """

        jobs = []
        for sample in self.samples:
            duplicate_file_prefix = os.path.join("alignment", sample.name, sample.name + ".sorted.dup.")
            input = duplicate_file_prefix + "bam"
            print_reads_output = duplicate_file_prefix + "recal.bam"
            base_recalibrator_output = duplicate_file_prefix + "recalibration_report.grp"

            jobs.append(concat_jobs([
                gatk.base_recalibrator(input, base_recalibrator_output),
                gatk.print_reads(input, print_reads_output, base_recalibrator_output),
                Job(input_files=[print_reads_output], output_files=[print_reads_output + ".md5"], command="md5sum " + print_reads_output + " > " + print_reads_output + ".md5")
            ], name="recalibration." + sample.name))

        report_file = os.path.join("report", "DnaSeq.recalibration.md")
        jobs.append(
            Job(
                [os.path.join("alignment", sample.name, sample.name + ".sorted.dup.recal.bam") for sample in self.samples],
                [report_file],
                command="""\
mkdir -p report && \\
cp \\
  {report_template_dir}/{basename_report_file} \\
  {report_file}""".format(
                    report_template_dir=self.report_template_dir,
                    basename_report_file=os.path.basename(report_file),
                    report_file=report_file
                ),
                report_files=[report_file],
                name="recalibration_report")
        )

        return jobs

########################################################################################################

    def metrics(self):
        """
        Compute metrics and generate coverage tracks per sample. Multiple metrics are computed at this stage:
        Number of raw reads, Number of filtered reads, Number of aligned reads, Number of duplicate reads,
        Median, mean and standard deviation of insert sizes of reads after alignment, percentage of bases
        covered at X reads (%_bases_above_50 means the % of exons bases which have at least 50 reads)
        whole genome or targeted percentage of bases covered at X reads (%_bases_above_50 means the % of exons
        bases which have at least 50 reads). A TDF (.tdf) coverage track is also generated at this step
        for easy visualization of coverage in the IGV browser.
        """

        ##check the library status
        library = {}
        for readset in self.readsets:
            if not library.has_key(readset.sample) :
                library[readset.sample]="SINGLE_END"
            if readset.run_type == "PAIRED_END" :
                library[readset.sample]="PAIRED_END"

        jobs = []
        for sample in self.samples:
            recal_file_prefix = os.path.join("alignment", sample.name, sample.name + ".sorted.dup.recal.")
            input = recal_file_prefix + "bam"

            job = picard.collect_multiple_metrics(input, recal_file_prefix + "all.metrics",  library_type=library[sample])
            job.name = "picard_collect_multiple_metrics." + sample.name
            jobs.append(job)

            # Compute genome coverage with GATK
            job = gatk.depth_of_coverage(input, recal_file_prefix + "all.coverage", bvatools.resolve_readset_coverage_bed(sample.readsets[0]))
            job.name = "gatk_depth_of_coverage.genome." + sample.name
            jobs.append(job)

            # Compute genome or target coverage with BVATools
            job = bvatools.depth_of_coverage(
                input,
                recal_file_prefix + "coverage.tsv",
                bvatools.resolve_readset_coverage_bed(sample.readsets[0]),
                other_options=config.param('bvatools_depth_of_coverage', 'other_options', required=False)
            )

            job.name = "bvatools_depth_of_coverage." + sample.name
            jobs.append(job)


            rm_igv=config.param('igvtools_compute_tdf', 'rm_igv_fol', required=False)

            if rm_igv.lower() in ['yes', 'true', 'y', '1']:
                job = concat_jobs([
                        igvtools.compute_tdf(input, input + ".tdf"),
                        Job(command="if [[ -d $HOME/igv ]]; then rm -r $HOME/igv; fi"),
                        Job(command="if [[ -f igv.log ]]; then rm igv.log; fi"),
                        ], name="igvtools_compute_tdf." + sample.name)
            else:
                job = igvtools.compute_tdf(input, input + ".tdf")
                job.name = "igvtools_compute_tdf." + sample.name

            jobs.append(job)

        return jobs

######################################################################################################################

    def picard_calculate_hs_metrics(self):
        """
        Compute on target percent of hybridisation based capture.
        """

        jobs = []

        created_interval_lists = []

        for sample in self.samples:
            coverage_bed = bvatools.resolve_readset_coverage_bed(sample.readsets[0])
            if coverage_bed:
                interval_list = re.sub("\.[^.]+$", ".interval_list", coverage_bed)

                if not interval_list in created_interval_lists:
                    job = tools.bed2interval_list(None, coverage_bed, interval_list)
                    job.name = "interval_list." + os.path.basename(coverage_bed)
                    jobs.append(job)
                    created_interval_lists.append(interval_list)

                recal_file_prefix = os.path.join("alignment", sample.name, sample.name + ".sorted.dup.recal.")
                job = picard.calculate_hs_metrics(recal_file_prefix + "bam", recal_file_prefix + "onTarget.tsv", interval_list)
                job.name = "picard_calculate_hs_metrics." + sample.name
                jobs.append(job)
        return jobs

#####################################################################################################################


    def gatk_callable_loci(self):
        """
        Computes the callable region or the genome as a bed track.
        """

        jobs = []

        for sample in self.samples:
            alignment_file_prefix = os.path.join("alignment", sample.name, sample.name + ".")

            job = gatk.callable_loci(alignment_file_prefix + "sorted.dup.recal.bam", alignment_file_prefix + "callable.bed", alignment_file_prefix + "callable.summary.txt")
            job.name = "gatk_callable_loci." + sample.name
            jobs.append(job)

        return jobs


#####################################################################################################################

    def extract_common_snp_freq(self):
        """
        Extracts allele frequencies of possible variants accross the genome.
        """

        jobs = []

        for sample in self.samples:
            alignment_file_prefix = os.path.join("alignment", sample.name, sample.name + ".")

            job = bvatools.basefreq(alignment_file_prefix + "sorted.dup.recal.bam", alignment_file_prefix + "commonSNPs.alleleFreq.csv", config.param('extract_common_snp_freq', 'common_snp_positions', type='filepath'), 0)
            job.name = "extract_common_snp_freq." + sample.name
            jobs.append(job)

        return jobs

##################################################################################################################

    def baf_plot(self):
        """
        Plots DepthRatio and B allele frequency of previously extracted alleles.
        """

        jobs = []

        for sample in self.samples:
            alignment_file_prefix = os.path.join("alignment", sample.name, sample.name + ".")

            job = bvatools.ratiobaf(alignment_file_prefix + "commonSNPs.alleleFreq.csv", alignment_file_prefix + "ratioBAF", config.param('baf_plot', 'common_snp_positions', type='filepath'))
            job.name = "baf_plot." + sample.name
            jobs.append(job)

        return jobs


##################################################################################################################

    def gatk_haplotype_caller(self):
        """
        GATK haplotype caller for snps and small indels.
        """

        jobs = []

        nb_haplotype_jobs = config.param('gatk_haplotype_caller', 'nb_jobs', type='posint')
        if nb_haplotype_jobs > 50:
            log.warning("Number of haplotype jobs is > 50. This is usually much. Anything beyond 20 can be problematic.")

        for sample in self.samples:
            alignment_directory = os.path.join("alignment", sample.name)
            haplotype_directory = os.path.join(alignment_directory, "rawHaplotypeCaller")
            input = os.path.join(alignment_directory, sample.name + ".sorted.dup.recal.bam")

            if nb_haplotype_jobs == 1:
                jobs.append(concat_jobs([
                    # Create output directory since it is not done by default by GATK tools
                    Job(command="mkdir -p " + haplotype_directory,removable_files=[haplotype_directory]),
                    gatk.haplotype_caller(input, os.path.join(haplotype_directory, sample.name + ".hc.g.vcf.bgz"))
                ], name="gatk_haplotype_caller." + sample.name))

            else:
                unique_sequences_per_job,unique_sequences_per_job_others = split_by_size(self.sequence_dictionary, nb_haplotype_jobs - 1)

                # Create one separate job for each of the first sequences
                for idx,sequences in enumerate(unique_sequences_per_job):
                    jobs.append(concat_jobs([
                        # Create output directory since it is not done by default by GATK tools
                        Job(command="mkdir -p " + haplotype_directory,removable_files=[haplotype_directory]),
                        gatk.haplotype_caller(input, os.path.join(haplotype_directory, sample.name + "." + str(idx) + ".hc.g.vcf.bgz"), intervals=sequences)
                    ], name="gatk_haplotype_caller." + sample.name + "." + str(idx)))

                # Create one last job to process the last remaining sequences and 'others' sequences
                jobs.append(concat_jobs([
                    # Create output directory since it is not done by default by GATK tools
                    Job(command="mkdir -p " + haplotype_directory,removable_files=[haplotype_directory]),
                    gatk.haplotype_caller(input, os.path.join(haplotype_directory, sample.name + ".others.hc.g.vcf.bgz"), exclude_intervals=unique_sequences_per_job_others)
                ], name="gatk_haplotype_caller." + sample.name + ".others"))

        return jobs

##################################################################################################################

    def merge_and_call_individual_gvcf(self):
        """
        Merges the gvcfs of haplotype caller and also generates a per sample vcf containing genotypes.
        """

        jobs = []
        nb_haplotype_jobs = config.param('gatk_haplotype_caller', 'nb_jobs', type='posint')

        for sample in self.samples:
            haplotype_file_prefix = os.path.join("alignment", sample.name, "rawHaplotypeCaller", sample.name)
            output_haplotype_file_prefix = os.path.join("alignment", sample.name, sample.name)
            if nb_haplotype_jobs == 1:
                gvcfs_to_merge = [haplotype_file_prefix + ".hc.g.vcf.bgz"]
            else:
                unique_sequences_per_job,unique_sequences_per_job_others = split_by_size(self.sequence_dictionary, nb_haplotype_jobs - 1)

                gvcfs_to_merge = [haplotype_file_prefix + "." + str(idx) + ".hc.g.vcf.bgz" for idx in xrange(len(unique_sequences_per_job))]
                gvcfs_to_merge.append(haplotype_file_prefix + ".others.hc.g.vcf.bgz")

            jobs.append(concat_jobs([
                gatk.cat_variants(gvcfs_to_merge, output_haplotype_file_prefix + ".hc.g.vcf.bgz"),
                gatk.genotype_gvcf([output_haplotype_file_prefix + ".hc.g.vcf.bgz"], output_haplotype_file_prefix + ".hc.vcf.bgz",config.param('gatk_merge_and_call_individual_gvcfs', 'options'))
            ], name="merge_and_call_individual_gvcf." + sample.name))

        return jobs

#################################################################################################################

    def combine_gvcf(self):
	"""
	Combine the per sample gvcfs of haplotype caller into one main file for all sample.
	"""
	jobs = []
	nb_haplotype_jobs = config.param('gatk_combine_gvcf', 'nb_haplotype', type='posint')
	nb_maxbatches_jobs = config.param('gatk_combine_gvcf', 'nb_batch', type='posint')

	# merge all sample in one shot
	if nb_maxbatches_jobs == 1 :
	    if nb_haplotype_jobs == 1:
		jobs.append(concat_jobs([
		    Job(command="mkdir -p variants"),
		    gatk.combine_gvcf([ os.path.join("alignment", sample.name, sample.name)+".hc.g.vcf.bgz" for sample in self.samples ], os.path.join("variants", "allSamples.hc.g.vcf.bgz"))],
		    name="gatk_combine_gvcf.AllSamples"))
	    else :
		unique_sequences_per_job,unique_sequences_per_job_others = split_by_size(self.sequence_dictionary, nb_haplotype_jobs - 1)

		# Create one separate job for each of the first sequences
		for idx,sequences in enumerate(unique_sequences_per_job):
		    obs.append(concat_jobs([
			Job(command="mkdir -p variants",removable_files=[os.path.join("variants", "allSamples") + "." + str(idx) + ".hc.g.vcf.bgz",os.path.join("variants", "allSamples") + "." + str(idx) + ".hc.g.vcf.bgz.tbi"]),
			gatk.combine_gvcf([ os.path.join("alignment", sample.name, sample.name)+".hc.g.vcf.bgz" for sample in self.samples ], os.path.join("variants", "allSamples") + "." + str(idx) + ".hc.g.vcf.bgz", intervals=sequences)
		    ], name="gatk_combine_gvcf.AllSample" + "." + str(idx)))

		# Create one last job to process the last remaining sequences and 'others' sequences
		job=gatk.combine_gvcf([ os.path.join("alignment", sample.name, sample.name)+".hc.g.vcf.bgz" for sample in self.samples ], os.path.join("alignment", "allSamples.others.hc.g.vcf.bgz"), exclude_intervals=unique_sequences_per_job_others)
		job.name="gatk_combine_gvcf.AllSample" + ".others"
		job.removable_files=[os.path.join("variants", "allSamples.others.hc.g.vcf.bgz"),os.path.join("variants", "allSamples.others.hc.g.vcf.bgz.tbi") ]
		jobs.append(job)
	else:
	    #Combine samples by batch (pre-defined batches number in ini)
	    sample_per_batch = int(math.ceil(len(self.samples)/float(nb_maxbatches_jobs)))
	    batch_of_sample = [ self.samples[i:(i+sample_per_batch)] for i in range(0,len(self.samples),sample_per_batch) ]
	    cpt = 0
	    batches = []
	    for batch in batch_of_sample :
		if nb_haplotype_jobs == 1:
		    jobs.append(concat_jobs([
			Job(command="mkdir -p variants",removable_files=[os.path.join("variants", "allSamples.batch" + str(cpt) + ".hc.g.vcf.bgz"),os.path.join("variants", "allSamples.batch" + str(cpt) + ".hc.g.vcf.bgz.tbi")]),
			gatk.combine_gvcf([ os.path.join("alignment", sample.name, sample.name)+".hc.g.vcf.bgz" for sample in batch ], os.path.join("variants", "allSamples.batch" + str(cpt) + ".hc.g.vcf.bgz"))
		    ], name="gatk_combine_gvcf.AllSamples.batch" + str(cpt)))
		else :
		    unique_sequences_per_job,unique_sequences_per_job_others = split_by_size(self.sequence_dictionary, nb_haplotype_jobs - 1)

		    # Create one separate job for each of the first sequences
		    for idx,sequences in enumerate(unique_sequences_per_job):
			jobs.append(concat_jobs([
			    Job(command="mkdir -p variants",removable_files=[os.path.join("variants", "allSamples") + ".batch" + str(cpt) + "." + str(idx) + ".hc.g.vcf.bgz",os.path.join("variants", "allSamples") + ".batch" + str(cpt) + "." + str(idx) + ".hc.g.vcf.bgz.tbi"]),
			    gatk.combine_gvcf([ os.path.join("alignment", sample.name, sample.name)+".hc.g.vcf.bgz" for sample in batch ], os.path.join("variants", "allSamples") + ".batch" + str(cpt) + "." + str(idx) + ".hc.g.vcf.bgz", intervals=sequences)
			], name="gatk_combine_gvcf.AllSample" + ".batch" + str(cpt) + "." + str(idx)))

		    # Create one last job to process the last remaining sequences and 'others' sequences
		    job=gatk.combine_gvcf([ os.path.join("alignment", sample.name, sample.name)+".hc.g.vcf.bgz" for sample in batch ], os.path.join("variants", "allSamples" + ".batch" + str(cpt) + ".others.hc.g.vcf.bgz"), exclude_intervals=unique_sequences_per_job_others)
		    job.name="gatk_combine_gvcf.AllSample" + ".batch" + str(cpt) + ".others"
		    job.removable_files=[os.path.join("variants", "allSamples" + ".batch" + str(cpt) + ".others.hc.g.vcf.bgz"),os.path.join("variants", "allSamples" + ".batch" + str(cpt) + ".others.hc.g.vcf.bgz.tbi")]
		    jobs.append(job)
		batches.append("batch" + str(cpt))
		cpt = cpt + 1

	    #Combine batches altogether
	    if nb_haplotype_jobs == 1:
		job=gatk.combine_gvcf([ os.path.join("variants", "allSamples." + batch_idx + ".hc.g.vcf.bgz") for batch_idx in batches ], os.path.join("variants", "allSamples.hc.g.vcf.bgz"))
		job.name="gatk_combine_gvcf.AllSamples.batches"
		jobs.append(job)
	    else :
		unique_sequences_per_job,unique_sequences_per_job_others = split_by_size(self.sequence_dictionary, nb_haplotype_jobs - 1)

		# Create one separate job for each of the first sequences
		for idx,sequences in enumerate(unique_sequences_per_job):
		    job=gatk.combine_gvcf([ os.path.join("variants", "allSamples." + batch_idx + "." + str(idx) + ".hc.g.vcf.bgz") for batch_idx in batches ], os.path.join("variants", "allSamples") + "." + str(idx) + ".hc.g.vcf.bgz", intervals=sequences)
		    job.name="gatk_combine_gvcf.AllSample" + "." + str(idx)
		    job.removable_files=[os.path.join("variants", "allSamples") + "." + str(idx) + ".hc.g.vcf.bgz",os.path.join("variants", "allSamples") + "." + str(idx) + ".hc.g.vcf.bgz.tbi"]
		    jobs.append(job)

		# Create one last job to process the last remaining sequences and 'others' sequences
		job=gatk.combine_gvcf([ os.path.join("variants", "allSamples." + batch_idx + ".others.hc.g.vcf.bgz") for batch_idx in batches ], os.path.join("variants", "allSamples" + ".others.hc.g.vcf.bgz"), exclude_intervals=unique_sequences_per_job_others)
		job.name="gatk_combine_gvcf.AllSample" + ".others"
		job.removable_files=[os.path.join("variants", "allSamples" + ".others.hc.g.vcf.bgz"),os.path.join("variants", "allSamples" + ".others.hc.g.vcf.bgz.tbi")]
		jobs.append(job)

        return jobs

####################################################################################################################

    def merge_and_call_combined_gvcf(self):
        """
        Merges the combined gvcfs and also generates a general vcf containing genotypes.
        """

        jobs = []
        nb_haplotype_jobs = config.param('gatk_combine_gvcf', 'nb_haplotype', type='posint')

        haplotype_file_prefix = os.path.join("variants","allSamples")
        output_haplotype = os.path.join("variants", "allSamples.hc.g.vcf.bgz")
        output_haplotype_genotyped = os.path.join("variants", "allSamples.hc.vcf.bgz")
        if nb_haplotype_jobs > 1:
            unique_sequences_per_job,unique_sequences_per_job_others = split_by_size(self.sequence_dictionary, nb_haplotype_jobs - 1)

            gvcfs_to_merge = [haplotype_file_prefix + "." + str(idx) + ".hc.g.vcf.bgz" for idx in xrange(len(unique_sequences_per_job))]
            gvcfs_to_merge.append(haplotype_file_prefix + ".others.hc.g.vcf.bgz")

            job = gatk.cat_variants(gvcfs_to_merge, output_haplotype)
            job.name = "merge_and_call_combined_gvcf.merge.AllSample"
            jobs.append(job)

        job = gatk.genotype_gvcf([output_haplotype], output_haplotype_genotyped ,config.param('gatk_merge_and_call_combined_gvcfs', 'options'))
        job.name = "merge_and_call_combined_gvcf.call.AllSample"
        jobs.append(job)

        return jobs

########################################################################################################

    def variant_recalibrator(self):
        """
        GATK VariantRecalibrator.
        The purpose of the variant recalibrator is to assign a well-calibrated probability to each variant call in a call set.
        You can then create highly accurate call sets by filtering based on this single estimate for the accuracy of each call.
        The approach taken by variant quality score recalibration is to develop a continuous, covarying estimate of the relationship
        between SNP call annotations (QD, MQ, HaplotypeScore, and ReadPosRankSum, for example) and the probability that a SNP
        is a true genetic variant versus a sequencing or data processing artifact. This model is determined adaptively based
        on "true sites" provided as input, typically HapMap 3 sites and those sites found to be polymorphic on the Omni 2.5M SNP
        chip array. This adaptive error model can then be applied to both known and novel variation discovered in the call set
        of interest to evaluate the probability that each call is real. The score that gets added to the INFO field of each variant
        is called the VQSLOD. It is the log odds ratio of being a true variant versus being false under the trained Gaussian mixture model.
        Using the tranche file generated by the previous step the ApplyRecalibration walker looks at each variant's VQSLOD value
        and decides which tranche it falls in. Variants in tranches that fall below the specified truth sensitivity filter level
        have their filter field annotated with its tranche level. This will result in a call set that simultaneously is filtered
        to the desired level but also has the information necessary to pull out more variants for a higher sensitivity but a
        slightly lower quality level.
        """

        jobs = []

        #generate the recalibration tranche files
        output_directory = "variants"
        recal_snps_other_options = config.param('variant_recalibrator', 'tranch_other_options_snps')
        recal_indels_other_options = config.param('variant_recalibrator', 'tranch_other_options_indels')
        variant_recal_snps_prefix = os.path.join(output_directory, "allSamples.hc.snps")
        variant_recal_indels_prefix = os.path.join(output_directory, "allSamples.hc.indels")

        jobs.append(concat_jobs([
            Job(command="mkdir -p " + output_directory),
            gatk.variant_recalibrator( [os.path.join(output_directory, "allSamples.hc.vcf.bgz")], recal_snps_other_options, variant_recal_snps_prefix + ".recal", variant_recal_snps_prefix + ".tranches", variant_recal_snps_prefix + ".R"),
            gatk.variant_recalibrator( [os.path.join(output_directory, "allSamples.hc.vcf.bgz")], recal_indels_other_options, variant_recal_indels_prefix + ".recal", variant_recal_indels_prefix + ".tranches", variant_recal_indels_prefix + ".R")
        ], name="variant_recalibrator.tranch.allSamples"))


        #apply the recalibration
        apply_snps_other_options = config.param('variant_recalibrator', 'apply_other_options_snps')
        apply_indels_other_options = config.param('variant_recalibrator', 'apply_other_options_indels')
        variant_apply_snps_prefix = os.path.join(output_directory, "allSamples.hc.snps")
        variant_apply_indels_prefix = os.path.join(output_directory, "allSamples.hc.indels")

        jobs.append(concat_jobs([
            Job(command="mkdir -p " + output_directory),
            gatk.apply_recalibration( os.path.join(output_directory, "allSamples.hc.vcf.bgz"), variant_apply_snps_prefix + ".recal", variant_apply_snps_prefix + ".tranches", apply_snps_other_options, variant_apply_snps_prefix + "_raw_indels.genotyped.vqsr.vcf.bgz"),
            gatk.apply_recalibration( variant_apply_snps_prefix + "_raw_indels.genotyped.vqsr.vcf.bgz", variant_apply_indels_prefix + ".recal", variant_apply_indels_prefix + ".tranches", apply_indels_other_options, os.path.join(output_directory, "allSamples.hc.vqsr.vcf"))
        ], name="variant_recalibrator.apply.allSamples"))

        return jobs


# The following two steps (multi_to_individual_VQSR and individual_VQSR_GVCF_intersect) were added so
# that following VQSR filration on the multi-sample VCF, the rest of the Jacek filters would be applied
# to single-sample VCF files
#####################################################################################################

    def multi_to_individual_VQSR(self):
        """
	Breakdown the variants/allSamples.hc.vqsr.vcf file into individual sample files using "bcftools view" method
        """

	jobs = []
	input = "variants/allSamples.hc.vqsr.vcf"

	for sample in self.samples:
		jobs.append(concat_jobs([
		    bcftools.multiToIndivVQSR2(input, sample.name, os.path.join(os.path.join("variants",sample.name), sample.name+".hc.vqsr.vcf"))
		    ], name="multi_to_individual_VQSR"+sample.name))

	return jobs



#####################################################################################################
    def individual_VQSR_GVCF_intersect(self):
	"""
	This function was added to handle multi-sample run of the pipeline. The previous "variant_recalibrator" step which creates VQSR
	tranches, should be run on multiple samples. Then the following steps should be performed:
	1) Trim the "sample.hc.g.vcf" file for each sample so that "bcftools isec" can work with it. This involves removing all lines with
	   only <NON_REF> in the ALT field and removing the ",<NON_REF>" phrase from other lines that have allele(s) in the ALT field
	   The resulting file is "sample.hc.g.trimmed.vcf"
	2) Create a "sample.PASSvqsr.vcf" subset of these individual vqsr files, which contain only variants that have "PASS"ed the vqsr filter
	3) Find the intersection of "sample.PASSvqsr.vcf" file with the "sample.hc.g.trimmed.vcf" gvcf file of each sample. The latter file
	   contains the original DP, MQ, etc. parameters.
	"""

            
	jobs=[]

	for sample in self.samples:
	    alignment_directory = os.path.join("alignment", sample.name)
	    variant_directory = os.path.join("variants",sample.name)
	    input = os.path.join(variant_directory, sample.name + ".hc.vqsr.vcf")
            bcftools_view_options = "-O z" 


            jobs.append(concat_jobs([
		Job(input_files=[os.path.join(alignment_directory, sample.name+".hc.g.vcf.bgz")], output_files=[os.path.join(alignment_directory, sample.name+".hc.g.vcf")], command="zcat "+ os.path.join(alignment_directory, sample.name+".hc.g.vcf.bgz") + " > "+ os.path.join(alignment_directory, sample.name+".hc.g.vcf")),
                Job(input_files=[os.path.join(alignment_directory, sample.name+".hc.g.vcf")], output_files=[os.path.join(alignment_directory,sample.name+".hc.g.trimmed.vcf")], command="cat " + os.path.join(alignment_directory, sample.name)+".hc.g.vcf | grep -v $'\t<NON_REF>' | sed 's/,<NON_REF>//' > " + os.path.join(alignment_directory, sample.name)+".hc.g.trimmed.vcf"),
		bcftools.view(os.path.join(alignment_directory, sample.name)+".hc.g.trimmed.vcf", os.path.join(alignment_directory, sample.name)+".hc.g.trimmed.vcf.gz", bcftools_view_options),
		tabix.tabix_index(os.path.join(alignment_directory, sample.name)+".hc.g.trimmed.vcf.gz", "vcf"),

                Job(input_files=[input], output_files=[os.path.join(variant_directory,sample.name+".hc.PASSvqsr.vcf")], command="cat " + input + " | grep -v $'\tVQSRTranche' > " + os.path.join(variant_directory, sample.name)+".hc.PASSvqsr.vcf"),
		bcftools.view(os.path.join(variant_directory, sample.name)+".hc.PASSvqsr.vcf", os.path.join(variant_directory, sample.name)+".hc.PASSvqsr.vcf.gz", bcftools_view_options),
		tabix.tabix_index(os.path.join(variant_directory, sample.name)+".hc.PASSvqsr.vcf.gz", "vcf"),

		bcftools.intersect(os.path.join(variant_directory, sample.name)+".hc.PASSvqsr.vcf.gz", sample.name, os.path.join(variant_directory, sample.name)+".0002.vcf")

	        ], name="individual_VQSR_GVCF_intersect"))


	return jobs

#####################################################################################################


    def varfilter(self):
        """
        Calling/Filtering varaints.
        If pileup was used, varfilter is called using SAMTools.
        If mpileup was used, varfilter is called using BCFTools.
        """

	variants_directory = os.path.join(self.output_dir, "variants")

        jobs = []

        for sample in self.samples:
            bcftools_view_options = "-O z"
            variant_directory = os.path.join("variants", sample.name)
            ref_fasta = config.param("bcftools_norm", "genome_fasta", type="filepath")

	    # "0002.vcf" file is generated by the bcftools intersect function. More details in "intersect" mothod of /bfx/bcftools.py
	    input_hc_vqsr=os.path.join(variant_directory, sample.name+".0002.vcf")
	    output_raw_vcf=os.path.join(variant_directory, sample.name+".raw.vcf")

            jobs.append(concat_jobs([
                        samtools.bcftools_varfilter(input_hc_vqsr, os.path.join(variant_directory, sample.name+".init.vcf")),

                        bcftools.view(os.path.join(variant_directory, sample.name+".init.vcf"), os.path.join(variant_directory, sample.name+".init.vcf.gz"), bcftools_view_options),
                        tabix.tabix_index(os.path.join(variant_directory, sample.name+".init.vcf.gz"), "vcf"),

                        bcftools.norm(os.path.join(variant_directory, sample.name+".init.vcf.gz"), os.path.join(variant_directory, sample.name+".walker.vcf"), "-m-both"),
                        bcftools.norm(os.path.join(variant_directory, sample.name+".walker.vcf"), os.path.join(variants_directory, sample.name+".flt.vcf"), "-f " + ref_fasta)
                        ], name="bcftools_varfilter"))

        return jobs



##################################################################################################################################


    def update_master_variants(self):
        """
        Store called variants in the alignment phase of the pipeline to a variant DB.
        The DB is a flat file - however you can opt for a real DBMS such PostGreSQL or MariaDB.
        """

        jobs = []

        out_dir = os.path.join(self.output_dir, "variants/db")
	variants_directory = os.path.join(self.output_dir, "variants")
        vcf_folder = os.path.join(out_dir, "vcf_folder")
        variants_fol = config.param("update_master_variants", "variants_fol")
        vardb_file = config.param("update_master_variants", "vardb_name")
        vcflist_file = config.param("update_master_variants", "vcflist_name")

        sampleNames = [sample.name for sample in self.samples]

        in_vcfs = [os.path.join(variants_directory, sample.name + ".flt.vcf") for sample in self.samples]

        path_vcf = os.path.join(self.output_dir, "variants")

        sampleNames = ",".join(sampleNames)

        jobs.append(concat_jobs([
            Job(command="mkdir -p " + vcf_folder),
            Job([os.path.join(variants_fol, vardb_file), os.path.join(variants_fol, vcflist_file)],
                [os.path.join(out_dir, vardb_file), os.path.join(out_dir, vcflist_file)],
                command="cp " + os.path.join(variants_fol, vardb_file) + " " + os.path.join(variants_fol, vcflist_file) + " " + out_dir),
            #forge_tools.update_master_variants(sampleNames, in_vcfs, path_vcf, out_dir, vcf_folder)
            forge_tools.update_master_variants(sampleNames, in_vcfs, path_vcf, variants_fol, vcf_folder)
        ], name="update_master_variants"))

        return jobs




    def convert2annovar(self):
        """
        Convert the filtered VCF files for each sample into an annovar file for the subsequent steps.
        """

        jobs = []


	variants_directory = os.path.join(self.output_dir, "variants")

        for sample in self.samples:

            vcf_file = os.path.join(variants_directory, sample.name + ".flt.vcf")

            output_prefix = os.path.join("annotation", sample.name)
            out_file = os.path.join(output_prefix, sample.name + ".annovarIn.txt")

            jobs.append(concat_jobs([
                Job(command="mkdir -p " + output_prefix),
                annovar.convert2annovar(vcf_file, out_file)
            ], name="convert2annovar." + sample.name))

        return jobs



    def annovar_annotation(self):
        """
        Annotate the filtered VCF files using [Annovar](http://annovar.openbioinformatics.org/en/latest/).
        This function uses Annovar to annotate: genes, dbsnp, 1000g, lib26, sift, polyphen, lrt, mutationtaster, gerp, phastcons.
        """

        jobs = []

        for sample in self.samples:
            output_prefix = os.path.join("annotation", sample.name, sample.name)
            annovar_file = os.path.join(output_prefix + ".annovarIn.txt")
            # Annovar appends ends to the out file depending on the call
            # Gene files' suffixes will be:  ".annovar_out.exonic_variant_function" and ".annovar_out.variant_function"
            # Table annovar file will be: ".annovar_out.hg19_multianno.txt"
            out_file = os.path.join(output_prefix + ".annovar_out")
            dbsnp_ver = config.param("annotate_genes", "dbsnp_version")

            protocol="1000g2015aug_all,snp"+ dbsnp_ver+",ljb26_all,exac03,clinvar_20150629,phastConsElements46way"
	    operation="f,f,f,f,f,r"

            jobs.append(concat_jobs([
                annovar.annotate_genes(annovar_file, out_file, self.get_option(sample.name, "SPLICING_THRESHOLD", "annotate_genes")),
                annovar.annotate_genes(annovar_file, out_file + ".splicing_ext", self.get_option(sample.name, "SPLICING_THRESHOLD_EXT", "annotate_genes")),

                annovar.table_annovar(annovar_file, out_file, protocol, operation)
            ], name="annovar_annotation." + sample.name))

        return jobs



    def combine_annovar_files(self):
        """
        Combine the annovar annotations.
        """

        jobs = []

	variants_directory = os.path.join(self.output_dir, "variants")

        for sample in self.samples:

            vcf_file = os.path.join(variants_directory, sample.name + ".flt.vcf")
            output_prefix = os.path.join("annotation", sample.name, sample.name + ".annovar_out")

            # Annovar output files
            exonic_vars = os.path.join(output_prefix + ".exonic_variant_function")
            var_function = os.path.join(output_prefix + ".variant_function")
            var_fn_extended = os.path.join(output_prefix + ".splicing_ext.variant_function")
            table_file = os.path.join(output_prefix + ".hg19_multianno.txt")
            out_file = os.path.join(output_prefix + ".combined.vcf")

            job = annovar.combine_annovar_files(vcf_file, exonic_vars, var_function, var_fn_extended, table_file, out_file)
            job.name = "combine_annovar_files." + sample.name
            jobs.append(job)

        return jobs


    def normalize_vcf(self):
        """
        Normalize the vcf using vt
        """

        jobs = []

        for sample in self.samples:

            vcf_file = os.path.join("annotation", sample.name, sample.name + ".annovar_out.combined.vcf")

            job = vt.decompose_and_normalize_mnps(vcf_file, sample.name + ".annovar_out.normalized.vcf.gz")
            job.name = "normalize_vcf." + sample.name
            jobs.append(job)

        return jobs


    def gemini_annotate(self):
        """
        Annotate using SnpEff or VEP
        """

        jobs = []

        # TODO: make path locations configurable
        # TODO: add support for snpeff
        for sample in self.samples:

            vcf_file = os.path.join("annotation", sample.name, sample.name + ".annovar_out.normalized.vcf.gz")
            # TODO: don't name with vep
            out_file = os.path.join("annotation", sample.name, sample.name + ".vep-annotated.vcf")

            job = vep.annotate(vcf_file, out_file)
            job.name = "vep." + sample.name
            jobs.append(job)

        return jobs


    def gemini_load_to_db(self):
        """
        Load each vcf into a db
        """

        jobs = []

        for sample in self.samples:

            # TODO: don't name with vep
            vcf_file = os.path.join("annotation", sample.name, sample.name + ".vep-annotated.vcf")
            db_name = os.path.join("annotation", sample.name, sample.name + ".gemini.db")

            job = gemini.gemini_annotations(vcf_file, db_name, os.path.join("annotation", sample.name))
            job.name = "gemini_load_to_db." + sample.name
            jobs.append(job)

        return jobs


    def prev_seen(self):
        """
        Annotating previously seen variants.
        """

        jobs = []

        for sample in self.samples:
            vcf_file = os.path.join("annotation", sample.name, sample.name + ".annovar_out.combined.vcf")
            ps_file = os.path.join("annotation", sample.name, sample.name + ".combined.ps.vcf")
            vardb_file = config.param("update_master_variants", "vardb_name")
            vardb = os.path.join(self.output_dir, "variants/db", vardb_file)

            job = forge_tools.prev_seen(vcf_file, ps_file, vardb, self.get_option(sample.name, "PREV_SEEN_DETAILS_THRESHOLD", "prev_seen"))
            job.name = "prev_seen." + sample.name
            jobs.append(job)

        return jobs



    def gene_mutation_counts(self):
        """
        Annotating details about GENE.
        """

        jobs = []

        for sample in self.samples:
            vcf_file = os.path.join("annotation", sample.name, sample.name + ".combined.ps.vcf")
            gmf_file = os.path.join("annotation", sample.name, sample.name + ".combined.ps.gmf.vcf")

            job = forge_tools.gene_mutation_counts(vcf_file, gmf_file)
            job.name = "gene_mutation_counts." + sample.name
            jobs.append(job)

        return jobs



    def evs(self):
        """
        Annotating with EXOME VARIANT SERVER.
        """

        jobs = []

        for sample in self.samples:
            vcf_file = os.path.join("annotation", sample.name, sample.name + ".combined.ps.gmf.vcf")
            evs_file = os.path.join("annotation", sample.name, sample.name + ".combined.ps.gmf.evs.vcf")

            job = forge_tools.evs(vcf_file, evs_file)
            job.name = "evs." + sample.name
            jobs.append(job)

        return jobs



    def omim(self):
        """
        Annotating OMIM.
        """

        jobs = []

        for sample in self.samples:
            vcf_file = os.path.join("annotation", sample.name, sample.name + ".combined.ps.gmf.evs.vcf")
            omim_file = os.path.join("annotation", sample.name, sample.name + ".combined.ps.gmf.evs.omim.vcf")

            job = forge_tools.omim(vcf_file, omim_file)
            job.name = "omim." + sample.name
            jobs.append(job)

        return jobs


    # Helper function - get the value of the option "option_name"
    # First try to get it from the sample in the sample info file
    # Next try to get the value from the "ALL" value in the sample info file
    # Finally, try to get the value from the main config file
    # All else fails, return None
    def get_option(self, sample_name, option_name, ini_section="filter"):

        sample_options = {}
        all_options = {}

        if sample_name in self.sample_info.keys():
            sample_options = self.sample_info[sample_name]
        if "ALL" in self.sample_info.keys():
            all_options = self.sample_info["ALL"]

        log.warning("OPTION: " + option_name)
        log.warning("sample_options: " + str(sample_options))
        log.warning("all_options: " + str(all_options))

        if option_name in sample_options.keys():
            return sample_options[option_name]
        elif option_name in all_options.keys():
            return all_options[option_name]
        elif config.param(ini_section, option_name, required=False):
            return config.param(ini_section, option_name, required=False)

        return None



    def filter(self):
        """
        Filter low quality variants.
        """

        jobs = []

        for sample in self.samples:

            # Set up the options
            minReadCount = self.get_option(sample.name, "MIN_READ_COUNT")
            minAltCount = self.get_option(sample.name, "MIN_ALT_COUNT")
            minSNVReadRatio = self.get_option(sample.name, "MIN_SNV_READ_RATIO")
            minIndelReadRatio = self.get_option(sample.name, "MIN_INDEL_READ_RATIO")
            minQ = self.get_option(sample.name, "MIN_QUALITY")
            minMapQ = self.get_option(sample.name, "MIN_MAPQ")
            prevSeenRemoveThreshold = self.get_option(sample.name, "REMOVE_PREV_SEEN")
            MAFThreshold = self.get_option(sample.name, "MAFThreshold")
            filterSSE = self.get_option(sample.name, "FILTER_SSE")
            filterRandoms = int(self.get_option(sample.name, "FILTER_RANDOMS"))
            cnvMaxP = self.get_option(sample.name, "CNV_MAX_P")
            maxCNVs = self.get_option(sample.name, "MAX_CNVS")

	    minQD = self.get_option(sample.name, "MIN_QD")
	    minReadPosRankSumIndel = self.get_option(sample.name, "MIN_READ_POS_RANK_SUM_INDEL")
	    maxFSIndel = self.get_option(sample.name, "MAX_FS_INDEL")
	    minReadPosRankSumSNP = self.get_option(sample.name, "MIN_READ_POS_RANK_SUM_SNP")
	    maxFSSNP = self.get_option(sample.name, "MAX_FS_SNP")
	    minMapQSNP = self.get_option(sample.name, "MIN_MAPQ_SNP")
	    minMQRankSumSNP = self.get_option(sample.name, "MIN_MQRankSum_SNP")

            options="""{minReadCount}{minAltCount}{minSNVReadRatio}{minIndelReadRatio}{minQ}{minMapQ}{prevSeenRemoveThreshold}{MAFThreshold}{filterSSE}{filterRandoms}{cnvMaxP}{maxCNVs}{minQD}{minReadPosRankSumIndel}{maxFSIndel}{minReadPosRankSumSNP}{maxFSSNP}{minMapQSNP}{minMQRankSumSNP}""".format(
                minReadCount=" \\\n  --minReadCount " + minReadCount if minReadCount else "",
                minAltCount=" \\\n  --minAltCount " + minAltCount if minAltCount else "",
                minSNVReadRatio=" \\\n  --minSNVReadRatio " + minSNVReadRatio if minSNVReadRatio else "",
                minIndelReadRatio=" \\\n  --minIndelReadRatio " + minIndelReadRatio if minIndelReadRatio else "",
                minQ=" \\\n  --minQ " + minQ if minQ else "",
                minMapQ=" \\\n  --minMapQ " + minMapQ if minMapQ else "",
                prevSeenRemoveThreshold=" \\\n  --prevSeenThreshold " + prevSeenRemoveThreshold if prevSeenRemoveThreshold else "",
                MAFThreshold=" \\\n  --maxMAF " + MAFThreshold if MAFThreshold else "",
                filterSSE=" \\\n  --filterSSE" if filterSSE else "",
                filterRandoms=" \\\n  --filterRandoms" if filterRandoms else "",
                cnvMaxP=" \\\n  --cnvMaxP " + cnvMaxP if cnvMaxP else "",
                maxCNVs=" \\\n  --maxCNVs " + maxCNVs if maxCNVs else "",

		minQD=" \\\n --minQD " + minQD if minQD else "",
		minReadPosRankSumIndel=" \\\n --minReadPosRankSumIndel " + minReadPosRankSumIndel if minReadPosRankSumIndel else "",
		maxFSIndel=" \\\n --maxFSIndel " + maxFSIndel if maxFSIndel else "",
		minReadPosRankSumSNP=" \\\n --minReadPosRankSumSNP " + minReadPosRankSumSNP if minReadPosRankSumSNP else "",
		maxFSSNP=" \\\n --maxFSSNP " + maxFSSNP if maxFSSNP else "",
		minMapQSNP=" \\\n --minMapQSNP " + minMapQSNP if minMapQSNP else "",
		minMQRankSumSNP=" \\\n --minMQRankSumSNP " + minMQRankSumSNP if minMQRankSumSNP else ""
            )

            include_intronic = self.get_option(sample.name, "INCLUDE_INTRONIC")
            keepOnlyVarTypes = self.get_option(sample.name, "keepOnlyVarTypes")
            keepOnlyVarTypes = keepOnlyVarTypes + ",intronic" if include_intronic else keepOnlyVarTypes

	    light_options = options + " \\\n  --keepOnlyVarTypes \"" + keepOnlyVarTypes + "\"" if keepOnlyVarTypes else ""

            # TODO: CHECK -- what if keepOnlyVarTypes isn't defined?
            if (keepOnlyVarTypes):
                keepWideVarTypes = "\"" + keepOnlyVarTypes + ",synonymous SNV,exonic-splicing,UTR3,UTR5,UTR5-UTR3\""
            else:
                keepWideVarTypes = "\"synonymous SNV,exonic-splicing,UTR3,UTR5,UTR5-UTR3\""

            exten_options = options + " \\\n  --keepOnlyVarTypes " + keepWideVarTypes


            # Set up the filter job
            output_prefix = os.path.join("annotation", sample.name, sample.name)
            vcf_file = os.path.join(output_prefix + ".combined.ps.gmf.evs.omim.vcf")
            light_output = os.path.join(output_prefix + ".combined.ps.gmf.evs.omim.light.flt.vcf")
            output = os.path.join(output_prefix + ".combined.ps.gmf.evs.omim.flt.vcf")
            all_output = os.path.join(output_prefix + ".combined.ps.gmf.evs.omim.all.flt.vcf")

            jobs.append(concat_jobs([
                forge_tools.filter(vcf_file, light_output, light_options),
                forge_tools.filter(vcf_file, output, exten_options),
                forge_tools.num_private_variants(sample.name, vcf_file, output_prefix + ".num_private_variants.txt"),
                forge_tools.filter(vcf_file, all_output, exten_options + " \\\n  --doNotRemoveFiltered")
            ], name="filter." + sample.name))

        return jobs



    def homozygosity(self):
        """
        Annotating homozygosity.
        """

        jobs = []

        for sample in self.samples:
            output_prefix = os.path.join("annotation", sample.name, sample.name)
            vcf_file = os.path.join(output_prefix + ".combined.ps.gmf.evs.omim.flt.vcf")
            light_vcf = os.path.join(output_prefix + ".combined.ps.gmf.evs.omim.light.flt.vcf")
            all_vcf = os.path.join(output_prefix + ".combined.ps.gmf.evs.omim.all.flt.vcf")

            output_hom = os.path.join(output_prefix + ".combined.flt.hom.vcf")
            light_hom = os.path.join(output_prefix + ".combined.light.flt.hom.vcf")
            all_hom = os.path.join(output_prefix + ".combined.all.flt.hom.vcf")

            roh_file = os.path.join(output_prefix + ".regions_of_homozygosity.txt")
            output_roh = os.path.join(output_prefix + ".annotatedROH.vcf")
            light_roh = os.path.join(output_prefix + ".annotatedROH.light.vcf")
            all_roh = os.path.join(output_prefix + ".annotatedROH.all.vcf")

            jobs.append(concat_jobs([
                forge_tools.homozygosity(vcf_file, output_hom),
                forge_tools.homozygosity(light_vcf, light_hom),
                forge_tools.homozygosity(all_vcf, all_hom),

                forge_tools.predict_roh(all_hom, roh_file),

                forge_tools.roh(output_hom, output_roh, roh_file),
                forge_tools.roh(light_hom, light_roh, roh_file),
                forge_tools.roh(all_hom, all_roh, roh_file)
            ], name="homozygosity." + sample.name))

        return jobs



    def allele_ratio_metrics(self):
        """
        Calculating ratio metrics and histogram.
        """

        jobs = []

        for sample in self.samples:
            output_prefix = os.path.join("annotation", sample.name, sample.name)
            vcf_file = os.path.join(output_prefix + ".annotatedROH.all.vcf")
            allele_metrics = os.path.join(output_prefix + ".allele_ratios.txt")
            allele_log = os.path.join(output_prefix + ".allele_ratios.log")

            jobs.append(concat_jobs([
                forge_tools.allele_ratio_metrics(vcf_file, allele_metrics),
                forge_tools.allele_ratio_hist(allele_metrics, output_prefix, allele_log)
            ], name="allele_ratio_metrics." + sample.name))

        return jobs



    def vcf2columns(self):
        """
        Converting annotated VCF to tab-separated values for Excel.
        """

        jobs = []

        # Modified to add the PP2_HDIV, PP2_HVAR, Exac, and Clinvar columns
        cols = "\"CHRPOS,VT,REF,ALT,ALTC,RDC,HMZ,ROH,QUAL,MQ,FILTER,PC,GENE,GMF,PSN,PS,DTLS,ID,THGMAF,EVSMAF,EVSGTC,EVSRD,EVSGRA,GERP,PHC,SIFT,PP2_HDIV,PP2_HVAR,MT,CADD_Phred,LRT,EXAC,CLINVAR,SNVAVG,GN,OMIM\""

        col_headers = "\"Position,,Variation,,Ref,,Alt,,#alt bases,,#reads,,Homozygosity,,In ROH,,VariantQ,,MapQ,,Filter,,Protein Change,,Gene,,Gene mutation frequency in our internal control samples (gene rank, #rare mutations, #truncating mutations),,#prev samples,,Prev seen in samples,,Info,,rsID,,MAF from 1000genomes,,EVS MAF,,EVS Genotype Counts (hom alt, het, hom ref),,EVS Avg Read Depth,,EVS Grantham,,GERP Score,,Phast cons score,,SIFT score (P.damaging),,Polyphen2_HDIV score (P.damaging),,Polyphen2_HVAR score (P.damaging),,MutationTaster score,,CADD_Phred,,LRT Score,,Exac frequency,,Clinvar,,NS SNV Avg Score,,Gene Description,,OMIM\""
        escape_cols = "\"GMF,EVSGTC\""

        for sample in self.samples:
            output_prefix = os.path.join("annotation", sample.name, sample.name)
            vcf_file = os.path.join(output_prefix + ".annotatedROH.vcf")
            light_vcf = os.path.join(output_prefix + ".annotatedROH.light.vcf")
            all_vcf = os.path.join(output_prefix + ".annotatedROH.all.vcf")

            output = os.path.join(output_prefix + ".combined_variants.txt")
            light_output = os.path.join(output_prefix + ".combined_variants.light.txt")
            all_output = os.path.join(output_prefix + ".combined_variants.all.txt")

            jobs.append(concat_jobs([
                forge_tools.vcf2columns(vcf_file, output, cols, escape_cols, col_headers),
                forge_tools.vcf2columns(light_vcf, light_output, cols, escape_cols, col_headers),
                forge_tools.vcf2columns(all_vcf, all_output, cols, escape_cols, col_headers)
            ], name="vcf2columns." + sample.name))

        return jobs



    def ucsc_link(self):
        """
        Annotate UCSC link.
        """

        jobs = []

        for sample in self.samples:
            output_prefix = os.path.join("annotation", sample.name, sample.name)
            std_input = os.path.join(output_prefix + ".combined_variants.txt")
            light_input = os.path.join(output_prefix + ".combined_variants.light.txt")
            all_input = os.path.join(output_prefix + ".combined_variants.all.txt")

            output = os.path.join(output_prefix + ".combined_variants.excel.txt")
            light_output = os.path.join(output_prefix + ".combined_variants.light.excel.txt")
            all_output = os.path.join(output_prefix + ".combined_variants.all.excel.txt")

            jobs.append(concat_jobs([
                forge_tools.ucsc_link(std_input, output),
                forge_tools.ucsc_link(light_input, light_output),
                forge_tools.ucsc_link(all_input, all_output)
            ], name="ucsc_link." + sample.name))

        return jobs


    @property
    def steps(self):
        return [
            self.picard_sam_to_fastq, #1
            # Part one of the pipeline [alignment_pipeline.pl]

            self.trimmomatic, #2
            self.merge_trimmomatic_stats, #3
            self.bwa_mem_picard_sort_sam, #4
            self.picard_merge_sam_files, #5
            self.gatk_indel_realigner, #6
            self.merge_realigned, #7
            self.fix_mate_by_coordinate, #8
            self.picard_mark_duplicates, #9
            self.recalibration, #10
            self.metrics, #11
            self.picard_calculate_hs_metrics, #12
            self.gatk_callable_loci, #13
            self.extract_common_snp_freq, #14
            self.baf_plot, #15
            self.gatk_haplotype_caller, #16
            self.merge_and_call_individual_gvcf, #17
            self.combine_gvcf, #18
            self.merge_and_call_combined_gvcf, #19
            self.variant_recalibrator, #20

            self.multi_to_individual_VQSR, #21
            self.individual_VQSR_GVCF_intersect, #22

            self.varfilter, #23

            # Part two of the pipeline
            self.update_master_variants, #24
            # Part three of the pipeline [annotation_pipeline.pl]
            self.convert2annovar, #25
            self.annovar_annotation, #26
            self.combine_annovar_files, #27
            self.prev_seen, #28
            self.gene_mutation_counts, #29
            self.evs, #30
            self.omim, #31
            self.filter, #32
            self.homozygosity, #33
            self.allele_ratio_metrics, #34
            self.vcf2columns, #35
            self.ucsc_link #36
        ]

if __name__ == '__main__':
    Forge()
