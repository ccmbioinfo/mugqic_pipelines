#!/usr/bin/env python


# MUGQIC Modules
from core.config import config
from core.job import Job


def flash(fastq1, fastq2):
    return Job(
        # TODO: specify input and output (when pipeline is more mature)
        # input_files=[fastq1, fastq2],
        # output_files=NotImplementedError,
        module_entries=['flash', 'module_flash'],
        command='flash {options} ' \
                '{fastq1} {fastq2}'.format(options=config.param('flash', 'options'),
                                           fastq1=fastq1,
                                           fastq2=fastq2))
