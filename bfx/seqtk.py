#!/usr/bin/env

from core.job import Job


def fastq_to_fasta(fastq, fasta, name='fastq_to_fasta'):
    return Job(name=name,
               input_files=[fastq],
               output_files=[fasta],
               module_entries=[['seqtk', 'module_seqtk']],
               command='seqtk seq -a {fastq} > {fasta}'.format(fastq=fastq,
                                                               fasta=fasta))
