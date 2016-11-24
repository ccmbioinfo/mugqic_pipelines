#!/usr/bin/env python
"""
Take only those reads from cow1 and cow2 with matching ids

Input:
remove_rrna/cow{1,2}_qual_unique_n_rRNA.fastq

Output:
remove_host_reads/cow{1,2}_matching_ids.fastq
"""

from Bio import SeqIO
import argparse


def normalize(id):
    """
    Eg. @SRR594215.10015644/1 ->
        @SRR594215.10015644

    Since the paired-end IDs are marked with a trailing '/1' or '/2',
        we need to remove this to compare IDs
    """
    return id.split('/')[0]


arg_parser = argparse.ArgumentParser()
arg_parser.add_argument('--fastq1')
arg_parser.add_argument('--fastq2')
arg_parser.add_argument('--output1')
arg_parser.add_argument('--output2')
args = arg_parser.parse_args()


reads1 = set(SeqIO.parse(args.fastq1, 'fastq'))
reads2 = set(SeqIO.parse(args.fastq2, 'fastq'))

ids_in_both = {normalize(read.id) for read in reads1} & {normalize(read.id) for read in reads2}

reads1 = {read for read in reads1 if normalize(read.id) in ids_in_both}
reads2 = {read for read in reads2 if normalize(read.id) in ids_in_both}

SeqIO.write(reads1, args.output1, 'fastq')
SeqIO.write(reads2, args.output2, 'fastq')
