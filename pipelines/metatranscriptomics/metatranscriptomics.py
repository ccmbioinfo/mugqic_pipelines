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

log = logging.getLogger(__name__)


# TODO: specify input/output files
class Metatranscriptomics(common.Illumina):
    # Location for pipeline scripts, to be used by all steps
    # 'metatranscriptomics/scripts'
    script_path = os.path.join(os.path.dirname(__file__), 'scripts')

    def format_fastq_headers(self):
        jobs = []

        output_dir = 'format_reads'

        for readset in self.readsets:
            output1 = join(output_dir, 'cow1_new.fastq')
            output2 = join(output_dir, 'cow2_new.fastq')

            jobs.append(Job(name='format_fastq_headers.' + readset.name,
                            module_entries=[['DEFAULT', 'module_perl']],
                            command='perl {script_path}/main_add_subID_reads_fastq.pl '
                                    '{input1} {output1} '
                                    '{input2} {output2}'.format(script_path=self.script_path,
                                                                input1=readset.fastq1,
                                                                output1=output1,
                                                                input2=readset.fastq2,
                                                                output2=output2)))

        return [Job(command='mkdir {}'.format(output_dir))].extend(jobs)

    def trimmomatic(self):
        input_dir = 'format_reads'
        input1 = join(input_dir, 'cow1_new.fastq')
        input2 = join(input_dir, 'cow2_new.fastq')

        output_dir = 'format_reads'
        output_paired1 = join(output_dir, 'cow1_qual_paired.fastq')
        output_unpaired1 = join(output_dir, 'cow1_qual_unpaired.fastq')
        output_paired2 = join(output_dir, 'cow2_qual_paired.fastq')
        output_unpaired2 = join(output_dir, 'cow2_qual_unpaired.fastq')

        job = trimmomatic.trimmomatic(input1,
                                      input2,
                                      output_paired1,
                                      output_unpaired1,
                                      output_paired2,
                                      output_unpaired2,
                                      None,
                                      None,
                                      adapter_file=config.param('trimmomatic', 'adapter_fasta'),
                                      trim_log=join(output_dir, 'cow.trim.log'))
        job.name = 'trimmomatic.cow'
        return job

    def flash(self):
        input = dict()
        input['dir'] = 'format_reads'
        input['fastq1'] = join(input['dir'], 'cow1_qual_paired.fastq')
        input['fastq2'] = join(input['dir'], 'cow2_qual_paired.fastq')

        output = dict()
        output['dir'] = 'format_reads'
        output['prefix'] = 'cow_qual'

        output['notCombined1'] = join(output['dir'], output['prefix'] + '.notCombined_1.fastq')
        output['notCombined2'] = join(output['dir'], output['prefix'] + '.notCombined_2.fastq')
        output['merged'] = join(output['dir'], output['prefix'] + '.extendedFrags.fastq')

        output[1] = join(output['dir'], 'cow1_qual_all.fastq')
        output[2] = join(output['dir'], 'cow2_qual_all.fastq')

        return [concat_jobs([flash.merge_short_reads(input['fastq1'], input['fastq2'], output['dir'], output['prefix']),
                             Job(command='cat {merged} {notCombined1} > {output1}'.format(merged=output['merged'],
                                                                                          notCombined1=output[
                                                                                              'notCombined1'],
                                                                                          output1=output[1])),
                             Job(command='cp {notCombined2} {output2}'.format(notCombined2=output['notCombined2'],
                                                                              output2=output[2]))],
                            name='flash')]

    def fastq_to_fasta(self):
        return [concat_jobs([Job(command='mkdir -p remove_duplicates'),
                             Job(input_files=['flash/cow1_qual_all.fastq',
                                              'flash/cow2_qual_all.fastq'],
                                 output_files=['remove_duplicates/cow1_qual_all.fasta',
                                               'remove_duplicates/cow2_qual_all.fasta'],
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
            self.flash,
            self.fastq_to_fasta,
            self.remove_duplicates,  # 5
            self.remove_abundant_rrna,
        ]


if __name__ == '__main__':
    Metatranscriptomics()
