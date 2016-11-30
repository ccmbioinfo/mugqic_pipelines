#!/usr/bin/env python


# MUGQIC Modules
from core.config import config
from core.job import Job


def merge_short_reads(fastq1, fastq2, output_dir, output_prefix):
    return Job(
        # TODO: specify input and output (when pipeline is more mature)
        # input_files=[fastq1, fastq2],
        # output_files=NotImplementedError,
        module_entries=['flash', 'module_flash'],
        command='flash -d {output_dir} -o {output_prefix} {options} ' \
                '{fastq1} {fastq2}'.format(output_dir=output_dir,
                                           output_prefix=output_prefix,
                                           options=config.param('flash', 'options'),
                                           fastq1=fastq1,
                                           fastq2=fastq2))
