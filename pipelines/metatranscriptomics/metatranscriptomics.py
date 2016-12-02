#!/usr/bin/env python

# Python Standard Modules
import logging
import os.path
from os.path import join
import sys

# Append mugqic_pipelines directory to Python library path
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0])))))

# MUGQIC Modules
from core.config import config
from core.job import Job, concat_jobs
from pipelines import common
from bfx import trimmomatic
from bfx import flash
from bfx import seqtk

log = logging.getLogger(__name__)


class Metatranscriptomics(common.Illumina):
    """
    Based off of the pipeline located here:
    http://bioinformatics-ca.github.io/analysis_of_metagenomic_data_module6_lab_2016/

    There are 4 phases to this pipeline so far:
    * format_reads
        Format read headers,
        trim reads,
        merge reads together
    * filter_reads
        Remove duplicates (temporarily),
        remove rRNA,
        remove reads from host
    * contig_assembly
        Assemble reads into larger contigs,
        map the reads to these contigs
    * search
        Search known microbial sequences using both contigs and singleton reads,
        try bwa, blat, and diamond
    * Further
    Each phase is associated with its own directory
    Within each phase directory, each sample/readset will have its own directory
    """

    # Location for pipeline scripts, to be used by all steps
    # 'metatranscriptomics/scripts'
    script_path = os.path.join(os.path.dirname(__file__), 'scripts')

    def format_fastq_headers(self):
        """
        Mark the headers of the fastq files with /1 or /2 to differentiate the paired-end reads

        Call 'main_add_subID_reads_fastq.pl'
        """
        jobs = []

        output_prefix = 'format_reads'
        for readset in self.readsets:
            output_dir = join(output_prefix, readset.name)

            output1 = join(output_dir, readset.name + '.1.formatted.fastq')
            output2 = join(output_dir, readset.name + '.2.formatted.fastq')

            jobs.append(concat_jobs([Job(command='mkdir {}'.format(output_dir)),
                                     Job(input_files=[readset.fastq1, readset.fastq2],
                                         output_files=[output1, output2],
                                         module_entries=[['DEFAULT', 'module_perl']],
                                         command='perl {script_path}/main_add_subID_reads_fastq.pl '
                                                 '{input1} {output1} '
                                                 '{input2} {output2}'.format(script_path=self.script_path,
                                                                             input1=readset.fastq1,
                                                                             output1=output1,
                                                                             input2=readset.fastq2,
                                                                             output2=output2))],
                                    name='format_fastq_headers.' + readset.name))

        return jobs

    def trimmomatic(self):
        jobs = []

        input_prefix = 'format_reads'
        output_prefix = 'format_reads'

        def get_inputs(readset):
            """
            :return: 2 fastq filenames for paired-end reads
            """
            input_dir = join(input_prefix, readset.name)
            return input_dir, \
                   join(input_dir, readset.name + '.1.formatted.fastq'), \
                   join(input_dir, readset.name + '.2.formatted.fastq')

        def get_outputs(readset):
            """
            :return: output directory name,
                     4 fastq filenames
            """
            output_dir = join(output_prefix, readset.name)
            return output_dir, \
                   join(output_dir, readset.name + '.1.qual_paired.fastq'), \
                   join(output_dir, readset.name + '.1.qual_unpaired.fastq'), \
                   join(output_dir, readset.name + '.2.qual_paired.fastq'), \
                   join(output_dir, readset.name + '.2.qual_unpaired.fastq')

        for readset in self.readsets:
            input_dir, input1, input2 = get_inputs(readset)
            output_dir, output_paired1, output_unpaired1, output_paired2, output_unpaired2 = get_outputs(readset)

            job = trimmomatic.trimmomatic(input1,
                                          input2,
                                          output_paired1,
                                          output_unpaired1,
                                          output_paired2,
                                          output_unpaired2,
                                          None,
                                          None,
                                          adapter_file=config.param('trimmomatic', 'adapter_fasta'),
                                          trim_log=join(output_prefix, readset.name + '.trim.log'))
            job.name = 'trimmomatic.' + readset.name
            jobs.append(job)

        return jobs

    def merge_overlapping_reads(self):
        """
        Reads from the paired-end fastqs are merged together.

        The merged reads are added to the first paired-end fastq
        """
        jobs = []

        input_prefix = 'format_reads'
        output_prefix = 'format_reads'

        def get_inputs(readset):
            """
            :return: 2 fastq filenames for paired-end reads
            """
            input_dir = join(input_prefix, readset.name)
            return join(input_dir, readset.name + '.1.qual_paired.fastq'), \
                   join(input_dir, readset.name + '.2.qual_paired.fastq')

        def get_flash_params(readset):
            """
            :return: flash's output directory and output prefix
            """
            return readset.name, join(output_prefix, readset.name)

        def get_flash_outputs(output_dir, flash_output_prefix):
            """
            Get the filenames that flash will output

            Flash will output 3 files, one for both fastqs, and one for the merged reads
            :return: 2 uncombined fastq files + merged reads fastq file
            """
            return join(output_dir, flash_output_prefix + '.notCombined_1.fastq'), \
                   join(output_dir, flash_output_prefix + '.notCombined_2.fastq'), \
                   join(output_dir, flash_output_prefix + '.extendedFrags.fastq')

        def get_outputs(readset):
            """
            :return: 2 output filenames
            """
            output_dir = join(output_prefix, readset.name)
            return join(output_dir, readset.name + '.1.qual_all.fastq'), \
                   join(output_dir, readset.name + '.2.qual_all.fastq')

        for readset in self.readsets:
            # Get the filenames
            input1, input2 = get_inputs(readset)
            flash_output_dir, flash_output_prefix = get_flash_params(readset)
            flash1, flash2, flash_merged = get_flash_outputs(flash_output_dir, flash_output_prefix)
            output1, output2 = get_outputs(readset)

            # Create jobs
            flash_job = flash.merge_overlapping_reads(input1, input2, flash_output_dir, flash_output_prefix)
            # Put the merged reads into fastq1
            fastq1_job = Job(input_files=[flash1, flash_merged],
                             output_files=[output1],
                             command='cat {flash1} {merged} > {output1}'.format(flash1=flash1, merged=flash_merged,
                                                                                output1=output1))
            # Rename fastq2 to be consistent with fastq1
            fastq2_job = Job(input_files=[flash2],
                             output_files=[output2],
                             command='mv {flash2} {output2}'.format(flash2=flash2, output2=output2))
            jobs.append(concat_jobs([flash_job, fastq1_job, fastq2_job], name='flash.' + readset.name))

        return jobs

    def fastq_to_fasta(self):
        jobs = []

        input_prefix = 'format_reads'
        output_prefix = 'filter_reads'

        for readset in self.readsets:
            input_dir = join(input_prefix, readset.name)
            fastq1 = join(input_dir, readset.name + '.1.qual_all.fastq')
            fastq2 = join(input_dir, readset.name + '.2.qual_all.fastq')

            output_dir = join(output_prefix, readset.name)
            fasta1 = join(output_dir, readset.name + '.1.qual_all.fasta')
            fasta2 = join(output_dir, readset.name + '.2.qual_all.fasta')

            jobs.append(concat_jobs([seqtk.fastq_to_fasta(fastq1, fasta1),
                                     seqtk.fastq_to_fasta(fastq2, fasta2)],
                                    name='fastq_to_fasta.' + readset.name))

        return jobs

    def remove_duplicates(self):
        return [
            Job(input_files=['remove_duplicates/cow1_qual_all.fasta',
                             'remove_duplicates/cow2_qual_all.fasta'],
                output_files=['remove_duplicates/cow1_qual_all_unique.fasta',
                              'remove_duplicates/cow1_qual_all_unique.uc',
                              'remove_duplicates/cow2_qual_all_unique.fasta',
                              'remove_duplicates/cow2_qual_all_unique.uc'],
                module_entries=[
                    ['remove_duplicates', 'module_usearch']
                ],
                command='''\
usearch --derep_fullseq --cluster remove_duplicates/cow1_qual_all.fasta --seedsout remove_duplicates/cow1_qual_all_unique.fasta --sizeout -uc remove_duplicates/cow1_qual_all_unique.uc
usearch --derep_fullseq --cluster remove_duplicates/cow2_qual_all.fasta --seedsout remove_duplicates/cow2_qual_all_unique.fasta --sizeout -uc remove_duplicates/cow2_qual_all_unique.uc
/hpf/largeprojects/ccmbio/nreinhardt/mugqic_pipelines/pipelines/metatranscriptomics/scripts/main_get_derepli_IDs.py''',
                name='remove_duplicates.cow')
        ]

    def remove_abundant_rrna(self):
        return [
            concat_jobs([
                Job(command='mkdir -p remove_abundant_rrna'),
                Job(input_files=['remove_duplicates/cow1_qual_all_unique.fasta'],
                    output_files=['remove_abundant_rrna/cow1_rRNA.log',
                                  'remove_abundant_rrna/cow1_rRNa.infernalout'],
                    module_entries=[
                        ['remove_abundant_rrna', 'module_infernal'],
                        ['remove_abundant_rrna', 'module_perl']
                    ],
                    command='cmscan -o remove_abundant_rrna/cow1_rRNA.log '
                            '--tblout remove_abundant_rrna/cow1_rRNA.infernalout '
                            '--noali --notextw --rfam -E 0.001 '
                            '{rfam} '
                            'remove_duplicates/cow1_qual_all_unique.fasta'.format(
                        rfam=config.param('remove_abundant_rrna', 'rfam_location'))),
                Job(input_files=['remove_duplicates/cow2_qual_all_unique.fasta'],
                    output_files=['remove_abundant_rrna/cow2_rRNA.log',
                                  'remove_abundant_rrna/cow2_rRNa.infernalout'],
                    module_entries=[
                        ['remove_abundant_rrna', 'module_infernal'],
                        ['remove_abundant_rrna', 'module_perl']
                    ],
                    command='cmscan -o remove_abundant_rrna/cow2_rRNA.log '
                            '--tblout remove_abundant_rrna/cow2_rRNA.infernalout '
                            '--noali --notextw --rfam -E 0.001 '
                            '{rfam} '
                            'remove_duplicates/cow2_qual_all_unique.fasta'.format(
                        rfam=config.param('remove_abundant_rrna', 'rfam_location')))
            ], name='remove_abundant_rrna.cow'
            )
        ]

    @property
    def steps(self):
        return [
            self.format_fastq_headers,
            self.trimmomatic,
            self.merge_overlapping_reads,
            self.fastq_to_fasta,
            self.remove_duplicates,  # 5
            self.remove_abundant_rrna,
        ]


if __name__ == '__main__':
    Metatranscriptomics()
