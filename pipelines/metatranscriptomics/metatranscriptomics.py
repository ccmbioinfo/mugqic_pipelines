#!/usr/bin/env python

# Python Standard Modules
import logging
import os
import sys

# Append mugqic_pipelines directory to Python library path
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0])))))

# MUGQIC Modules
from core.config import config
from core.job import Job, concat_jobs
from pipelines import common
from bfx import trimmomatic


log = logging.getLogger(__name__)


class Metatranscriptomics(common.Illumina):
    """
    Metatranscriptomics
    ===============================
    """
    def __init__(self):
        # Add pipeline specific arguments
        # self._lastPGStep = {}
        super(Metatranscriptomics, self).__init__()

    def format_fastq(self):
        return [Job(name='format_fastq.cow',
                    input_files=['../pipeline-testing/cow1.fastq',
                                 '../pipeline-testing/cow2.fastq'],  # TODO: use select_input_files, candidate_files
                    output_files=['../pipeline-testing/cow1_new.fastq',
                                  '../pipeline-testing/cow2_new.fastq'],
                    module_entries=[['DEFAULT', 'module_perl']],
                    command='perl {script_path}/scripts/main_add_subID_reads_fastq.pl'
                            '../pipeline-testing/cow'.format(script_path=config.param('DEFAULT', 'metatranscriptomics_location')))]

    def trimmomatic(self):
        return [concat_jobs([Job(command='mkdir -p ' + 'trim'),
                             trimmomatic.trimmomatic('../pipeline-testing/cow1_new.fastq',
                                                     '../pipeline-testing/cow2_new.fastq',
                                                     'trim/cow1_qual_paired.fastq',
                                                     'trim/cow1_qual_unpaired.fastq',
                                                     'trim/cow2_qual_paired.fastq',
                                                     'trim/cow2_qual_unpaired.fastq',
                                                     None,
                                                     None,
                                                     adapter_file=config.param('trimmomatic', 'adapter_fasta'),
                                                     trim_log='trim/cow.trim.log')],
                            name='trimmomatic.cow')]

    def flash(self):
        return [concat_jobs([Job(command='mkdir -p ' + 'flash'),
                             Job(input_files=['trim/cow1_qual_paired.fastq',
                                              'trim/cow2_qual_paired.fastq'],
                                 output_files=['flash/cow_qual.extendedFrags.fastq',
                                               'flash/cow_qual.hist',
                                               'flash/cow_qual.histogram',
                                               'flash/cow_qual.notCombined_1.fastq',
                                               'flash/cow_qual.notCombined_2.fastq'], # TODO
                                 command='{flash} -M 75 -p 64 -t 2 -o cow_qual -d flash '
                                         'trim/cow1_qual_paired.fastq trim/cow2_qual_paired.fastq'.format(flash=config.param('flash', 'location'))),
                             Job(input_files=['flash/cow_qual.extendedFrags.fastq',
                                              'flash/cow_qual.notCombined_1.fastq'],
                                 output_files=['flash/cow1_qual_all.fastq'],
                                 command='cat flash/cow_qual.extendedFrags.fastq flash/cow_qual.notCombined_1.fastq > flash/cow1_qual_all.fastq'),
                             Job(input_files=['flash/cow_qual.notCombined_2.fastq'],
                                 output_files=['flash/cow2_qual_all.fastq'],
                                 command='cp flash/cow_qual.notCombined_2.fastq flash/cow2_qual_all.fastq')],
                            name='flash.cow')]

    def fastq_to_fasta(self):
        return [concat_jobs([Job(command='mkdir -p remove_duplicates'),
                             Job(input_files=['flash/cow1_qual_all.fastq',
                                              'flash/cow2_qual_all.fastq'],
                                 output_files=['fastq_to_fasta'],
                                 module_entries=[
                                     ['fastq_to_fasta', 'module_seqtk'],
                                     ['fastq_to_fasta', 'module_perl']
                                 ],
                                 command='''\
seqtk seq -a flash/cow1_qual_all.fastq > remove_duplicates/cow1_qual_all.fasta
seqtk seq -a flash/cow2_qual_all.fastq > remove_duplicates/cow2_qual_all.fasta''')],
                            name='fastq_to_fasta.cow')]

    def remove_duplicates(self):
        return [
            Job(input_files=['remove_duplicates/cow1_qual_all.fastq',
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
usearch --derep_fullseq --cluster remove_duplicates/cow2_qual_all.fasta --seedsout remove_duplicates/cow2_qual_all_unique.fasta --sizeout -uc remove_duplicates/cow2_qual_all_unique.uc''',
                name='remove_duplicates.cow')
        ]

    @property
    def steps(self):
        return [
            self.format_fastq,
            self.trimmomatic,
            self.flash,
            self.fastq_to_fasta,
            self.remove_duplicates
        ]

if __name__ == '__main__':
    Metatranscriptomics()
