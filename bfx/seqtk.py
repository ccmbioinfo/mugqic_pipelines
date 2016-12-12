#!/usr/bin/env

from core.job import Job, concat_jobs, mkdir


def fastq_to_fasta(fastq, fasta, job_name='fastq_to_fasta'):
    return concat_jobs([mkdir(fasta),
                        Job(input_files=[fastq],
                            output_files=[fasta],
                            module_entries=[['seqtk', 'module_seqtk']],
                            command='seqtk seq -a {fastq} > {fasta}'.format(fastq=fastq,
                                                                            fasta=fasta))],
                       name=job_name)
