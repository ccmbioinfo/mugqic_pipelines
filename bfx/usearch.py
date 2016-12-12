#!/usr/bin/env python

from core.config import config
from core.job import Job


def cluster_duplicates(input_fasta, output_fasta, output_uc, job_name='cluster_duplicates'):
    return Job(name=job_name,
               input_files=[input_fasta],
               output_files=[output_fasta, output_uc],
               module_entries=[['DEFAULT', 'module_usearch']],
               command='usearch --derep_fullseq --cluster {input_fasta} '
                       '--seedsout {output_fasta} --sizeout '
                       '--uc {output_uc}'.format(input_fasta=input_fasta,
                                                 output_fasta=output_fasta,
                                                 output_uc=output_uc))