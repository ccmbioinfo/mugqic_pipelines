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
from bfx import picard
from bfx import forge_tools
from bfx import bvatools
from bfx import samtools
from bfx import bcftools
from bfx import tabix
from bfx import varscan
from bfx import annovar

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


    def fastqc1(self):
        """
        Produces read quality metrics using fastqc
        """
        jobs = []

        for readset in self.readsets:
            alignment_directory = os.path.join("alignment", readset.sample.name)
            fastqc_directory = os.path.join(alignment_directory, "fastqc", readset.name)

            # Find input readset FASTQs from original FASTQs in the readset sheet
            if readset.run_type == "PAIRED_END":
                candidate_input_files = [[readset.fastq1, readset.fastq2]]
                if readset.bam:
                    candidate_input_files.append([re.sub("\.bam$", ".pair1.fastq.gz", readset.bam), re.sub("\.bam$", ".pair2.fastq.gz", readset.bam)])
                fastq_files = self.select_input_files(candidate_input_files)
            elif readset.run_type == "SINGLE_END":
                candidate_input_files = [[readset.fastq1]]
                if readset.bam:
                    candidate_input_files.append([re.sub("\.bam$", ".single.fastq.gz", readset.bam)])
                fastq_files = self.select_input_files(candidate_input_files)
            else:
                raise Exception("Error: run type \"" + readset.run_type + 
                "\" is invalid for readset \"" + readset.name + "\" (should be PAIRED_END or SINGLE_END)!")

            jobs.append(concat_jobs([
                Job(command="mkdir -p " + fastqc_directory),
                fast.fastqc(fastq_files, fastqc_directory) 
            ], name="fastqc1." + readset.name))

        return jobs

    
    def fastqc2(self):
        """
        Produces read quality metrics of already TRIMMED reads using fastqc. 
        """
        
        jobs = []
        for readset in self.readsets:
            trim_file_prefix = os.path.join("trim", readset.sample.name, readset.name + ".trim.")
            alignment_directory = os.path.join("alignment", readset.sample.name)
            fastqc_directory = os.path.join(alignment_directory, "fastqc", readset.name)

            # Find input readset FASTQs from previous trimmomatic job
            # Does not check previous step fastqs because this step is made specifically to check trimmed fastqs
            if readset.run_type == "PAIRED_END":
                fastq_files = [trim_file_prefix + "pair1.fastq.gz", trim_file_prefix + "pair2.fastq.gz"]
            elif readset.run_type == "SINGLE_END":
                fastq_files = [trim_file_prefix + "single.fastq.gz"]
            else:
                raise Exception("Error: run type \"" + readset.run_type +
                "\" is invalid for readset \"" + readset.name + "\" (should be PAIRED_END or SINGLE_END)!")

            jobs.append(concat_jobs([
                Job(command="mkdir -p " + fastqc_directory),
                fast.fastqc(fastq_files, fastqc_directory)
            ], name="fastqc2." + readset.name))

        return jobs


    def bwa_aln(self):
        """
        Find the SA coordinates of the input reads using [BWA](http://bio-bwa.sourceforge.net/): bwa aln.
        This is done per sequencing readset.

        This step takes as input files:

        1. Trimmed FASTQ files if available
        2. Else, FASTQ files from the readset file if available
        3. Else, FASTQ output files from previous picard_sam_to_fastq conversion of BAM files
        """
        jobs = []

        for readset in self.readsets:
            pg_file = os.path.join("trim", readset.sample.name, readset.sample.name + ".pg.txt")
            trim_file_prefix = os.path.join("trim", readset.sample.name, readset.name + ".trim.")
            alignment_directory = os.path.join("alignment", readset.sample.name)
            aln_sai_prefix = os.path.join(alignment_directory, readset.name, readset.name)

            # Find input readset FASTQs first from previous trimmomatic job, then from original FASTQs in the readset sheet
            if readset.run_type == "PAIRED_END":
                candidate_input_files = [[trim_file_prefix + "pair1.fastq.gz", trim_file_prefix + "pair2.fastq.gz"]]
                if readset.fastq1 and readset.fastq2:
                    candidate_input_files.append([readset.fastq1, readset.fastq2])
                if readset.bam:
                    candidate_input_files.append([re.sub("\.bam$", ".pair1.fastq.gz", readset.bam), re.sub("\.bam$", ".pair2.fastq.gz", readset.bam)])
                [fastq1, fastq2] = self.select_input_files(candidate_input_files)
                aln_sai_prefix = aln_sai_prefix + ".pair1"
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
            
            # Do bwa aln on the two paired end files separately
            # Add this job if single end or paired end
            bwa_job = bwa.aln(fastq1, aln_sai_prefix + ".aln.sai")
            jobs.append(concat_jobs([
                Job(command="mkdir -p " + os.path.dirname(aln_sai_prefix)),
                bwa_job,
                forge_tools.add_pg_file(pg_file, "bwa_aln.1-" + readset.name, self.get_ver("bwa"), "\"" + bwa_job.command + "\"", self._lastPGStep[readset.name])
            ], name="bwa_aln.1." + readset.name))

            # Only add this job if paired end
            if fastq2:
                aln_sai2 = os.path.join(alignment_directory, readset.name, readset.name + ".pair2.aln.sai")
                bwa_job = bwa.aln(fastq2, aln_sai2)
                jobs.append(concat_jobs([
                    Job(command="mkdir -p " + os.path.dirname(aln_sai_prefix)),
                    bwa_job,
                    forge_tools.add_pg_file(pg_file, "bwa_aln.2-" + readset.name, self.get_ver("bwa"), "\"" + bwa_job.command + "\"", self._lastPGStep[readset.name])
                ], name="bwa_aln.2." + readset.name))

            self._lastPGStep[readset.name] = "bwa_aln.1-" + readset.name

        return jobs
##### COPIED REPORT SECTION FROM DNASEQ BWA_MEM. MAYBE WANT TO ADD A REPORT LATER
#         report_file = os.path.join("report", "DnaSeq.bwa_mem_picard_sort_sam.md")
#         jobs.append(
#             Job(
#                 [os.path.join("alignment", readset.sample.name, readset.name, readset.name + ".sorted.bam") for readset in self.readsets],
#                 [report_file],
#                 [["bwa_mem_picard_sort_sam", "module_pandoc"]],
#                 command="""\
# mkdir -p report && \\
# pandoc --to=markdown \\
#   --template {report_template_dir}/{basename_report_file} \\
#   --variable scientific_name="{scientific_name}" \\
#   --variable assembly="{assembly}" \\
#   {report_template_dir}/{basename_report_file} \\
#   > {report_file}""".format(
#                     scientific_name=config.param("bwa_mem_picard_sort_sam", "scientific_name"),
#                     assembly=config.param("bwa_mem_picard_sort_sam", "assembly"),
#                     report_template_dir=self.report_template_dir,
#                     basename_report_file=os.path.basename(report_file),
#                     report_file=report_file
#                 ),
#                 report_files=[report_file],
#                 name="bwa_mem_picard_sort_sam_report")
#         )


    def bwa_sam(self):
        """
        Generate alignments in the SAM format using [BWA](http://bio-bwa.sourceforge.net/): bwa samse/sampe
        This is done per sequencing readset.
        """
        jobs = []

        for readset in self.readsets:
            trim_file_prefix = os.path.join("trim", readset.sample.name, readset.name + ".trim.")
            alignment_directory = os.path.join("alignment", readset.sample.name)
            aln_sai_prefix = os.path.join(alignment_directory, readset.name, readset.name)
            readset_sam = os.path.join(alignment_directory, readset.name, readset.name + ".aln.sam")
            pg_file = os.path.join("trim", readset.sample.name, readset.sample.name + ".pg.txt")
            
            # Find input readset FASTQs first from previous trimmomatic job, then from original FASTQs in the readset sheet
            if readset.run_type == "PAIRED_END":
                candidate_input_files = [[trim_file_prefix + "pair1.fastq.gz", trim_file_prefix + "pair2.fastq.gz"]]
                if readset.fastq1 and readset.fastq2:
                    candidate_input_files.append([readset.fastq1, readset.fastq2])
                if readset.bam:
                    candidate_input_files.append([re.sub("\.bam$", ".pair1.fastq.gz", readset.bam), re.sub("\.bam$", ".pair2.fastq.gz", readset.bam)])
                [fastq1, fastq2] = self.select_input_files(candidate_input_files)

                # Get the command without the read_group parameter because the PG header line can't have another "@" tag in it
                pg_com = bwa.sampe(fastq1, fastq2, aln_sai_prefix + ".pair1.aln.sai", aln_sai_prefix + ".pair2.aln.sai", readset_sam).command

                # Call bwa sampe if paired end
                jobs.append(concat_jobs([
                    Job(command="mkdir -p " + os.path.dirname(readset_sam)),
                    bwa.sampe(fastq1, 
                              fastq2, 
                              aln_sai_prefix + ".pair1.aln.sai",
                              aln_sai_prefix + ".pair2.aln.sai",
                              readset_sam,
                              read_group="'@RG" + \
                                "\tID:" + readset.name + \
                                "\tSM:" + readset.sample.name + \
                                ("\tLB:" + readset.library if readset.library else "") + \
                                "\tPL:Illumina" + \
                                "'"
                    ),
                    forge_tools.add_pg_file(pg_file, "bwa_sampe-" + readset.name, self.get_ver("bwa"), "\"" + pg_com + "\"", self._lastPGStep[readset.name])
                ], name="bwa_sampe." + readset.name))

                self._lastPGStep[readset.name] = "bwa_sampe-" + readset.name

            elif readset.run_type == "SINGLE_END":
                candidate_input_files = [[trim_file_prefix + "single.fastq.gz"]]
                if readset.fastq1:
                    candidate_input_files.append([readset.fastq1])
                if readset.bam:
                    candidate_input_files.append([re.sub("\.bam$", ".single.fastq.gz", readset.bam)])
                [fastq1] = self.select_input_files(candidate_input_files)
                fastq2 = None
                
                # Get the command without the read_group parameters because the PG header line can't have another "@" tag in it
                pg_com = bwa.samse(fastq1, aln_sai_prefix + ".aln.sai", readset_sam)

                # Call bwa same if single end
                jobs.append(concat_jobs([
                    Job(command="mkdir -p " + os.path.dirname(readset_sam)),
                    bwa.samse(fastq1,
                              aln_sai_prefix + ".aln.sai",
                              readset_sam,
                              read_group="'@RG" + \
                                "\tID:" + readset.name + \
                                "\tSM:" + readset.sample.name + \
                                ("\tLB:" + readset.library if readset.library else "") + \
                                "\tPL:Illumina" + \
                                "'"
                    ),
                    forge_tools.add_pg_file(pg_file, "bwa_samse-" + readset.name, self.get_ver("bwa"), "\"" + pg_com + "\"", self._lastPGStep[readset.name])
                ], name="bwa_samse." + readset.name))

                self._lastPGStep[readset.name] = "bwa_samse-" + readset.name
            else:
                raise Exception("Error: run type \"" + readset.run_type +
                "\" is invalid for readset \"" + readset.name + "\" (should be PAIRED_END or SINGLE_END)!")

        return jobs


    def picard_sort_sam(self):
        """
        Convert SAM files to BAM and sorting using [Picard](http://broadinstitute.github.io/picard/).
        """
        jobs = []

        for readset in self.readsets:
            alignment_directory = os.path.join("alignment", readset.sample.name)
            readset_sam = os.path.join(alignment_directory, readset.name, readset.name + ".aln.sam")
            readset_bam = os.path.join(alignment_directory, readset.name, readset.name + ".sorted.bam")
            pg_file = os.path.join("trim", readset.sample.name, readset.sample.name + ".pg.txt")

            job = picard.sort_sam(readset_sam, readset_bam, "coordinate")
            jobs.append(concat_jobs([
                        job,
                        forge_tools.add_pg_file(pg_file, "picard_sort_sam-" + readset.name, self.get_ver("picard"), "\"" + job.command + "\"", self._lastPGStep[readset.name])
            ], name="picard_sort_sam." + readset.name))

            self._lastPGStep[readset.sample.name] = "picard_sort_sam-" + readset.name

        return jobs


    def picard_merge_sam_files(self):
        """
        Merging SAM files and producing single sorted BAM using [Picard](http://broadinstitute.github.io/picard/).

        Picard doesn't seem to work right when you try to merge SAM files into BAM and sort at the same time!
        BAM files end up much smaller than they should be with fewer reads.
        So we sort and produce BAMs first, then merge.
        """
        jobs = []

        for sample in self.samples:
            alignment_directory = os.path.join("alignment", sample.name)
            readset_bams = self.select_input_files([[os.path.join(alignment_directory, readset.name, readset.name + ".sorted.bam") for readset in sample.readsets], [readset.bam for readset in sample.readsets]])
            sample_bam = os.path.join(alignment_directory, sample.name + ".sorted.bam")
            pg_file = os.path.join("trim", sample.name, sample.name + ".pg.txt")

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
                picard_job = picard.merge_sam_files(readset_bams, sample_bam),

                job = concat_jobs([
                    mkdir_job,
                    picard_job,
                    forge_tools.add_pg_file(pg_file, "picard_merge_sam_files", self.get_ver("picard"), "\"" + picard_job.command + "\"", self._lastPGStep[sample.name])
                ], name="picard_merge_sam_files." + sample.name)

                self._lastPGStep[sample.name] = "picard_merge_sam_files"

            jobs.append(job)

        return jobs


    def gatk_indel_realigner(self):
        """
        Insertion and deletion realignment is performed on regions where multiple base mismatches
        are preferred over indels by the aligner since it can appear to be less costly by the algorithm.
        Such regions will introduce false positive variant calls which may be filtered out by realigning
        those regions properly. Realignment is done using [GATK](https://www.broadinstitute.org/gatk/).
        The reference genome is divided by a number regions given by the `nb_jobs` parameter.
        """
        jobs = []

        nb_jobs = config.param("gatk_indel_realigner", "nb_jobs", type="posint")
        if nb_jobs > 50:
            log.warning("Number of realign jobs is > 50. This is usually much. Anything beyond 20 can be problematic.")

        for sample in self.samples:
            alignment_directory = os.path.join("alignment", sample.name)
            realign_directory = os.path.join(alignment_directory, "realign")
            input = os.path.join(alignment_directory, sample.name + ".sorted.bam")
            pg_file = os.path.join("trim", sample.name, sample.name + ".pg.txt")

            if nb_jobs == 1:
                realign_prefix = os.path.join(realign_directory, "all")
                realign_intervals = realign_prefix + ".intervals"
                output_bam = realign_prefix + ".bam"
                sample_output_bam = os.path.join(alignment_directory, sample.name + ".realigned.bam")

                # Variables for PG File
                create_target_job = gatk.realigner_target_creator(input, realign_intervals)
                temp_lastPGStep = "gatk_realigner_target_creator"
                realigner_job = gatk.indel_realigner(input, output_bam, target_intervals=realign_intervals)

                jobs.append(concat_jobs([
                    Job(command="mkdir -p " + realign_directory, removable_files=[realign_directory]),
                    create_target_job,
                    realigner_job,
                    # Create sample realign symlink since no merging is required
                    Job([output_bam], [sample_output_bam], command="ln -s -f " + os.path.relpath(output_bam, os.path.dirname(sample_output_bam)) + " " + sample_output_bam),
                    forge_tools.add_pg_file(pg_file, "gatk_realigner_target_creator", self.get_ver("gatk"), "\"" + create_target_job.command + "\"", self._lastPGStep[sample.name]),
                    forge_tools.add_pg_file(pg_file, "gatk_indel_realigner", self.get_ver("gatk"), "\"" + realigner_job.command + "\"", temp_lastPGStep)
                ], name="gatk_indel_realigner." + sample.name))

                self._lastPGStep[sample.name] = "gatk_indel_realigner"

            else:
                # The first sequences are the longest to process.
                # Each of them must be processed in a separate job.
                unique_sequences_per_job = [sequence["name"] for sequence in self.sequence_dictionary[0:min(nb_jobs - 1, len(self.sequence_dictionary))]]

                # Create one separate job for each of the first sequences
                for sequence in unique_sequences_per_job:
                    realign_prefix = os.path.join(realign_directory, sequence)
                    realign_intervals = realign_prefix + ".intervals"
                    intervals=[sequence]
                    if unique_sequences_per_job.index(sequence) == 0:
                        intervals.append("unmapped")
                    output_bam = realign_prefix + ".bam"

                    # Variables for PG File
                    create_target_job = gatk.realigner_target_creator(input, realign_intervals, intervals=[sequence])
                    temp_lastPGStep = "gatk_realigner_target_creator-" + sequence
                    realigner_job = gatk.indel_realigner(input, output_bam, target_intervals=realign_intervals, intervals=intervals)

                    jobs.append(concat_jobs([
                        # Create output directory since it is not done by default by GATK tools
                        Job(command="mkdir -p " + realign_directory),
                        create_target_job,
                        realigner_job,
                        forge_tools.add_pg_file(pg_file, "gatk_realigner_target_creator-" + sequence, self.get_ver("gatk"), "\"" + create_target_job.command + "\"", self._lastPGStep[sample.name]),
                        forge_tools.add_pg_file(pg_file, "gatk_indel_realigner-" + sequence, self.get_ver("gatk"), "\"" + realigner_job.command + "\"", temp_lastPGStep)
                    ], name="gatk_indel_realigner." + sample.name + "." + sequence))

                # Create one last job to process the last remaining sequences and 'others' sequences
                realign_prefix = os.path.join(realign_directory, "others")
                realign_intervals = realign_prefix + ".intervals"
                output_bam = realign_prefix + ".bam"

                # Variables for PG File
                create_target_job = gatk.realigner_target_creator(input, realign_intervals, exclude_intervals=unique_sequences_per_job)
                temp_lastPGStep = "gatk_realigner_target_creator-others"
                realigner_job = gatk.indel_realigner(input, output_bam, target_intervals=realign_intervals, exclude_intervals=unique_sequences_per_job)

                jobs.append(concat_jobs([
                    # Create output directory since it is not done by default by GATK tools
                    Job(command="mkdir -p " + realign_directory, removable_files=[realign_directory]),
                    create_target_job,
                    realigner_job,
                    forge_tools.add_pg_file(pg_file, "gatk_realigner_target_creator-others", self.get_ver("gatk"), "\"" + create_target_job.command + "\"", self._lastPGStep[sample.name]),
                    forge_tools.add_pg_file(pg_file, "gatk_indel_realigner-others", self.get_ver("gatk"), "\"" + realigner_job.command + "\"", temp_lastPGStep)
                ], name="gatk_indel_realigner." + sample.name + ".others"))

                self._lastPGStep[sample.name] = "gatk_indel_realigner-others"

        return jobs


    def merge_realigned(self):
        """
        BAM files of regions of realigned reads are merged per sample using [Picard](http://broadinstitute.github.io/picard/).
        """
        jobs = []

        nb_jobs = config.param("gatk_indel_realigner", "nb_jobs", type="posint")

        for sample in self.samples:
            alignment_directory = os.path.join("alignment", sample.name)
            realign_directory = os.path.join(alignment_directory, "realign")
            merged_realigned_bam = os.path.join(alignment_directory, sample.name + ".realigned.bam")
            pg_file = os.path.join("trim", sample.name, sample.name + ".pg.txt")

            # if nb_jobs == 1, symlink has been created in indel_realigner and merging is not necessary
            if nb_jobs > 1:
                realigned_bams = [os.path.join(realign_directory, sequence["name"] + ".bam") for sequence in self.sequence_dictionary[0:min(nb_jobs - 1, len(self.sequence_dictionary))]]
                realigned_bams.append(os.path.join(realign_directory, "others.bam"))

                job = picard.merge_sam_files(realigned_bams, merged_realigned_bam)
                jobs.append(concat_jobs([
                            job,
                            forge_tools.add_pg_file(pg_file, "merge_realigned", self.get_ver("picard"), "\"" + job.command + "\"", self._lastPGStep[sample.name])
                ], name="merge_realigned." + sample.name))

                self._lastPGStep[sample.name] = "merge_realigned"

        report_file = os.path.join("report", "DnaSeq.gatk_indel_realigner.md")
        jobs.append(
            Job(
                [os.path.join("alignment", sample.name, sample.name + ".realigned.bam") for sample in self.samples],
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


    def recalibration(self):
        """
        Recalibrate base quality scores of sequencing-by-synthesis reads in an aligned BAM file. After recalibration, 
        the quality scores in the QUAL field in each read in the output BAM are more accurate in that
        the reported quality score is closer to its actual probability of mismatching the reference genome.
        Moreover, the reaclibration tool attempts to correct for variation in quality with machine cycle
        and sequence context, and by doing so, provides not only more accurate quality scores but also
        more widely dispersed ones.
        """
        jobs = []

        for sample in self.samples:
            alignment_directory = os.path.join("alignment", sample.name)
            merged_realigned_bam = os.path.join(alignment_directory, sample.name + ".realigned.bam")
            recal_table = os.path.join(alignment_directory, sample.name + ".recalibration.table")
            recalibrated_bam = os.path.join(alignment_directory, sample.name + ".recalibrated.bam")
            pg_file = os.path.join("trim", sample.name, sample.name + ".pg.txt")

            base_recal_job = gatk.base_recalibrator(merged_realigned_bam, recal_table)
            print_reads_job = gatk.print_reads(merged_realigned_bam, recalibrated_bam, recal_table)

            jobs.append(concat_jobs([
                        samtools.index(merged_realigned_bam),
                        base_recal_job,
                        print_reads_job,
                        Job([recalibrated_bam], [recalibrated_bam + ".md5"], command="md5sum " + recalibrated_bam + " > " + recalibrated_bam + ".md5"),
                        forge_tools.add_pg_file(pg_file, "gatk_base_recalibrator", self.get_ver("gatk"), "\"" + base_recal_job.command + "\"", self._lastPGStep[sample.name]),
                        forge_tools.add_pg_file(pg_file, "gatk_print_reads", self.get_ver("gatk"), "\"" + print_reads_job.command + "\"", "gatk_base_recalibrator")
            ], name="recalibration." + sample.name))

            self._lastPGStep[sample.name] = "gatk_print_reads"
            
        return jobs


    def picard_sort_recalibrated(self):
        """
        Sort the recalibrated (and merged) BAM files using [Picard](http://broadinstitute.github.io/picard/).
        """
        jobs = []

        for sample in self.samples:
            alignment_directory = os.path.join("alignment", sample.name)
            recalibrated_bam = os.path.join(alignment_directory, sample.name + ".recalibrated.bam")
            recalibrated_sorted_bam = os.path.join(alignment_directory, sample.name + ".recalibrated.sorted.bam")
            pg_file = os.path.join("trim", sample.name, sample.name + ".pg.txt")

            job = picard.sort_sam(recalibrated_bam, recalibrated_sorted_bam, "coordinate")
            jobs.append(concat_jobs([
                        job,
                        forge_tools.add_pg_file(pg_file, "picard_sort_recalibrated", self.get_ver("picard"), "\"" + job.command + "\"", self._lastPGStep[sample.name])
            ], name="picard_sort_recalibrated." + sample.name))

            self._lastPGStep[sample.name] = "picard_sort_recalibrated"

        return jobs


    def picard_mark_duplicates(self):
        """
        Mark duplicates. Aligned reads per sample are duplicates if they have the same 5' alignment positions
        (for both mates in the case of paired-end reads). All but the best pair (based on alignment score)
        will be marked as a duplicate in the BAM file. Marking duplicates is done using [Picard](http://broadinstitute.github.io/picard/).
        """
        jobs = []

        for sample in self.samples:
            alignment_file_prefix = os.path.join("alignment", sample.name, sample.name + ".")
            input = alignment_file_prefix + "recalibrated.sorted.bam"
            output = alignment_file_prefix + "recalibrated.sorted.dup.bam"
            metrics_file = alignment_file_prefix + "recalibrated.sorted.dup.metrics"
            pg_file = os.path.join("trim", sample.name, sample.name + ".pg.txt")

            job = picard.mark_duplicates([input], output, metrics_file)
            jobs.append(concat_jobs([
                        job,
                        forge_tools.add_pg_file(pg_file, "picard_mark_duplicates", self.get_ver("picard"), "\"" + job.command + "\"", self._lastPGStep[sample.name])
            ], name="picard_mark_duplicates." + sample.name))

            self._lastPGStep[sample.name] = "picard_mark_duplicates"

        report_file = os.path.join("report", "DnaSeq.picard_mark_duplicates.md")
        jobs.append(
            Job(
                [os.path.join("alignment", sample.name, sample.name + ".recalibrated.sorted.dup.bam") for sample in self.samples],
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


    def picard_fix_mate(self):
        """
        Fix the insert sizes and mate pair flags using [Picard](http://broadinstitute.github.io/picard/).
        """
        jobs = []

        for sample in self.samples:
            alignment_file_prefix = os.path.join("alignment", sample.name, sample.name + ".")
            input = alignment_file_prefix + "recalibrated.sorted.dup.bam"
            output = alignment_file_prefix + "recalibrated.sorted.dup.fixmate.bam"
            pg_file = os.path.join("trim", sample.name, sample.name + ".pg.txt")

            job = picard.fix_mate_information(input, output)
            jobs.append(concat_jobs([
                        job,
                        forge_tools.add_pg_file(pg_file, "picard_fix_mate", self.get_ver("picard"), "\"" + job.command + "\"", self._lastPGStep[sample.name])
            ], name="picard_fix_mate." + sample.name))

            self._lastPGStep[sample.name] = "picard_fix_mate"

        return jobs


    def update_header(self):
        """
        Update BAM file header to conform to FORGE specifications
        """
        jobs = []

        for sample in self.samples:
            alignment_file_prefix = os.path.join("alignment", sample.name, sample.name + ".")
            input = alignment_file_prefix + "recalibrated.sorted.dup.fixmate.bam"
            output =alignment_file_prefix + "recalibrated.sorted.dup.fixmate.final.bam" 

            pg_file = os.path.join("trim", sample.name, sample.name + ".pg.txt")
            header_file = alignment_file_prefix + "header.sam"

            job = forge_tools.update_bam_header(input, output, pg_file, header_file)
            job.name = "update_header." + sample.name
            jobs.append(job)

        return jobs


    def picard_collect_multiple_metrics(self):
        """ 
        Compute multiple metrics per sample using [Picard](http://broadinstitute.github.io/picard/).
        """

        # Check the library status
        library = {}
        for readset in self.readsets:
            if not library.has_key(readset.sample) :
                library[readset.sample]="SINGLE_END"
            if readset.run_type == "PAIRED_END" :
                library[readset.sample]="PAIRED_END"

        jobs = []
        
        for sample in self.samples:
            input = os.path.join("alignment", sample.name, sample.name + ".recalibrated.sorted.dup.fixmate.final.bam")
            final_metrics = os.path.join("alignment", sample.name, sample.name + ".dup.metrics")

            job = picard.collect_multiple_metrics(input, final_metrics, library_type=library[sample])
            job.name = "picard_collect_multiple_metrics." + sample.name
            jobs.append(job)

        return jobs

    def gatk_depth_of_coverage(self):
        """
        Compute the genome coverage using [GATK](https://www.broadinstitute.org/gatk/).
        """
        jobs = []

        for sample in self.samples:
            input = os.path.join("alignment", sample.name, sample.name + ".recalibrated.sorted.dup.fixmate.final.bam")
            coverage_directory = os.path.join("alignment", sample.name, "coverage")
            out_file_prefix = os.path.join(coverage_directory, sample.name + ".final.dup.genelist_coverage")

            jobs.append(concat_jobs([
                # Create output directory since it is not done by default by GATK tools
                Job(command="mkdir -p " + coverage_directory),
                gatk.depth_of_coverage(input, out_file_prefix, bvatools.resolve_readset_coverage_bed(sample.readsets[0])),
                forge_tools.aggregate_coverage(out_file_prefix + ".sample_interval_summary", out_file_prefix + ".txt"),
                forge_tools.aggregate_coverage(out_file_prefix + ".sample_interval_summary", out_file_prefix + "_with_intervals.txt", "--printIntervals")
            ], name="gatk_depth_of_coverage.genome." + sample.name))


            # Run on the BAM WITHOUT markdup
            coverage_dup_directory = os.path.join("alignment", sample.name, "coverage_including_dups")
            input = os.path.join("alignment", sample.name, sample.name + ".recalibrated.sorted.bam")
            out_file_prefix = os.path.join(coverage_dup_directory, sample.name + ".final.dup.genelist_coverage")
            
            jobs.append(concat_jobs([
                # Create output directory since it is not done by default by GATK tools
                Job(command="mkdir -p " + coverage_dup_directory),
                gatk.depth_of_coverage(input, out_file_prefix, bvatools.resolve_readset_coverage_bed(sample.readsets[0])),
                forge_tools.aggregate_coverage(out_file_prefix + ".sample_interval_summary", out_file_prefix + ".txt"),
                forge_tools.aggregate_coverage(out_file_prefix + ".sample_interval_summary", out_file_prefix + "_with_intervals.txt", "--printIntervals")
             ], name="gatk_depth_of_coverage.genome_with_dup." + sample.name))

        # # In the mcgill pipeline - all the intermediate BAM files are unlinked
        # unlink "${noArchiveSampleFolder}/${sampleName}.t30l30.sorted.bam";
        # unlink "${noArchiveSampleFolder}/${sampleName}.t30l30.sorted.bai";
        # unlink "${noArchiveSampleFolder}/${sampleName}.t30l30.realigned.bam";
        # unlink "${noArchiveSampleFolder}/${sampleName}.t30l30.realigned.bai";
        # unlink "${noArchiveSampleFolder}/${sampleName}.t30l30.realigned.sorted.bam";
        # unlink "${noArchiveSampleFolder}/${sampleName}.t30l30.realigned.sorted.bai";
        # unlink "${noArchiveSampleFolder}/${sampleName}.t30l30.realigned.markdup.sorted.bam";
        # unlink "${noArchiveSampleFolder}/${sampleName}.t30l30.realigned.markdup.sorted.bai";
        # unlink "${noArchiveSampleFolder}/${sampleName}.t30l30.realigned.markdup.sorted.fixmate.bam";
        # unlink "${noArchiveSampleFolder}/${sampleName}.t30l30.realigned.markdup.sorted.fixmate.bai";
        # unlink glob "${noArchiveSampleFolder}/*.sai";
        # unlink glob "${noArchiveSampleFolder}/*.fastq.gz";
            
        return jobs


    def generate_approximate_windows(self, nb_jobs):
        if nb_jobs <= len(self.sequence_dictionary):
            return [sequence["name"] + ":1-" + str(sequence["length"]) for sequence in self.sequence_dictionary]
        else:
            total_length = sum([sequence["length"] for sequence in self.sequence_dictionary])
            approximate_window_size = int(math.floor(total_length / (nb_jobs - len(self.sequence_dictionary))))
            windows = []

            for sequence in self.sequence_dictionary:
                for start, end in [[pos, min(pos + approximate_window_size - 1, sequence["length"])] for pos in range(1, sequence["length"] + 1, approximate_window_size)]:
                    windows.append(sequence["name"] + ":" + str(start) + "-" + str(end))
            return windows


    def mpileup(self):
        """
        Running SAMTOOLS mpileup. One mpileup file created per sample/chromosome
        """
        jobs = []

        for sample in self.samples:
            bcftools_call_options = "-vmO v"
            input_bam = os.path.join("alignment", sample.name, sample.name + ".recalibrated.sorted.dup.fixmate.final.bam")
            output = os.path.join("variants", sample.name + ".raw.vcf")


            mpileup_job = pipe_jobs([
                    samtools.mpileup([input_bam], None, config.param("samtools_mpileup", "mpileup_other_options"), regionFile=bvatools.resolve_readset_coverage_bed(sample.readsets[0])),
                    samtools.bcftools_call("-", None, bcftools_call_options),
                    ])

            
            jobs.append(concat_jobs([
                        Job(command="mkdir -p variants"),
                        pipe_jobs([
                                mpileup_job,
                                forge_tools.add_to_vcf_header(input_bam, output, mpileup_job.command)
                                ])], name="samtools_mpileup." + sample.name))

        return jobs


    def varfilter(self):
        """
        Calling/Filtering varaints.
        If pileup was used, varfilter is called using SAMTools.
        If mpileup was used, varfilter is called using BCFTools.
        """
        
        jobs = []

        for sample in self.samples:
            bcftools_view_options = "-O z"
            variant_prefix = os.path.join("variants", sample.name)
            ref_fasta = config.param("bcftools_norm", "genome_fasta", type="filepath")

            jobs.append(concat_jobs([
                        samtools.bcftools_varfilter(variant_prefix + ".raw.vcf", variant_prefix + ".init.vcf"),
                        bcftools.view(variant_prefix + ".init.vcf", variant_prefix + ".init.vcf.gz", bcftools_view_options),
                        tabix.tabix_index(variant_prefix + ".init.vcf.gz", "vcf"),
                        bcftools.norm(variant_prefix + ".init.vcf.gz", variant_prefix + ".walker.vcf", "-m-both"),
                        bcftools.norm(variant_prefix + ".walker.vcf", variant_prefix + ".flt.vcf", "-f " + ref_fasta)
                        ], name="bcftools_varfilter"))

            return jobs


    def varscan(self):
        """
        Calling variants using Varscan.
        """

        jobs = []
        
        for sample in self.samples:
            input_bam = os.path.join("alignment", sample.name, sample.name + ".recalibrated.sorted.dup.fixmate.final.bam")        
            variant_file_prefix = os.path.join("variants", sample.name)
            mpileup_file = variant_file_prefix + ".mpileup"
            vcf_prefix = os.path.join("variants", sample.name + ".flt.vcf")

            jobs.append(concat_jobs([
                samtools.mpileup([input_bam], mpileup_file, config.param("varscan", "mpileup_other_options")),
                varscan.mpileupcns_jacek(mpileup_file, variant_file_prefix + ".varscan.out", config.param("varscan", "other_options", required=False))
            ], name="varscan." + sample.name))

        return jobs

            
    def update_master_variants(self):
        """
        Store called variants in the alignment phase of the pipeline to a variant DB.
        The DB is a flat file - however you can opt for a real DBMS such PostGreSQL or MariaDB.
        """
        
        jobs = []

        out_dir = os.path.join(self.output_dir, "variants/db")
        vcf_folder = os.path.join(out_dir, "vcf_folder")
        variants_fol = config.param("update_master_variants", "variants_fol")
        vardb_file = config.param("update_master_variants", "vardb_name")
        vcflist_file = config.param("update_master_variants", "vcflist_name")

        sampleNames = [sample.name for sample in self.samples]
#        in_vcfs = [os.path.join(self.output_dir, "variants", sample.name + ".flt.vcf") for sample in self.samples]
        in_vcfs = [os.path.join("variants", sample.name + ".flt.vcf") for sample in self.samples]
        path_vcf = os.path.join(self.output_dir, "variants")
        sampleNames = ",".join(sampleNames)

        jobs.append(concat_jobs([
            Job(command="mkdir -p " + vcf_folder),
            Job([os.path.join(variants_fol, vardb_file), os.path.join(variants_fol, vcflist_file)], 
                [os.path.join(out_dir, vardb_file), os.path.join(out_dir, vcflist_file)], 
                command="cp " + os.path.join(variants_fol, vardb_file) + " " + os.path.join(variants_fol, vcflist_file) + " " + out_dir),
            forge_tools.update_master_variants(sampleNames, in_vcfs, path_vcf, out_dir, vcf_folder)
        ], name="update_master_variants"))

        return jobs


    def convert2annovar(self):
        """
        Convert the filtered VCF files for each sample into an annovar file for the subsequent steps.
        """
 
        jobs = []

        for sample in self.samples:
            vcf_file = os.path.join("variants", sample.name + ".flt.vcf")
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
            dbsnp_ver = config.param("DEFAULT", "dbsnp_version")
            protocol="1000g2015aug_all,snp"+ dbsnp_ver+",ljb26_all,phastConsElements46way"
            operation="f,f,f,r"

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
        
        for sample in self.samples:
            vcf_file = os.path.join("variants", sample.name + ".flt.vcf")
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
        
        # Run the filtering on each sample
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

            options="""{minReadCount}{minAltCount}{minSNVReadRatio}{minIndelReadRatio}{minQ}{minMapQ}{prevSeenRemoveThreshold}{MAFThreshold}{filterSSE}{filterRandoms}{cnvMaxP}{maxCNVs}""".format(
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
                maxCNVs=" \\\n  --maxCNVs " + maxCNVs if maxCNVs else ""
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

        # Maybe modify the columns and column headers because now have PP2_HDIV + PP2_HVAR
        cols = "\"CHRPOS,VT,REF,ALT,ALTC,RDC,HMZ,ROH,QUAL,DP4,PV4,MQ,FILTER,PC,GENE,GMF,PSN,PS,DTLS,ID,THGMAF,EVSMAF,EVSGTC,EVSRD,EVSGRA,GERP,PHC,SIFT,PP2,MT,LRT,SNVAVG,GN,OMIM\""
        col_headers = "\"Position,,Variation,,Ref,,Alt,,#alt bases,,#reads,,Homozygosity,,In ROH,,VariantQ,,DP4 (fwd ref, rev ref, fwd alt, rev alt),,PV4 (strandBias,baseQBias,mapQBias,endDistBias) p-values,,MapQ,,Filter,,Protein Change,,Gene,,Gene mutation frequency in our internal control samples (gene rank, #rare mutations, #truncating mutations),,#prev samples,,Prev seen in samples,,Info,,rsID,,MAF from 1000genomes,,EVS MAF,,EVS Genotype Counts (hom alt, het, hom ref),,EVS Avg Read Depth,,EVS Grantham,,GERP Score,,Phast cons score,,SIFT score (P.damaging),,Polyphen2 score (P.damaging),,MutationTaster score,,LRT Score,,NS SNV Avg Score,,Gene Description,,OMIM\""
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
            self.picard_sam_to_fastq,
            # Part one of the pipeline [alignment_pipeline.pl]
            self.fastqc1,
            self.trimmomatic,
            self.merge_trimmomatic_stats,
            self.fastqc2,
            self.bwa_aln,
            self.bwa_sam,
            self.picard_sort_sam,
            self.picard_merge_sam_files,
            self.gatk_indel_realigner,
            self.merge_realigned,
            self.recalibration,
            self.picard_sort_recalibrated,
            self.picard_mark_duplicates,
            self.picard_fix_mate,
            self.update_header,
            self.picard_collect_multiple_metrics,
            self.gatk_depth_of_coverage,
            self.mpileup,
            self.varfilter,
            self.varscan,
            # Part two of the pipeline
            self.update_master_variants,
            # Part three of the pipeline [annotation_pipeline.pl]
            self.convert2annovar,
            # gene, dbsnp, 1000g, lib26, sift, polyphen, lrt, mutationtaster, gerp, phastcons
            self.annovar_annotation, 
            self.combine_annovar_files,
            self.prev_seen,
            self.gene_mutation_counts,
            self.evs,
            self.omim,
            self.filter,
            self.homozygosity,
            self.allele_ratio_metrics,
            self.vcf2columns,
            self.ucsc_link
        ]

if __name__ == '__main__':
    Forge()
