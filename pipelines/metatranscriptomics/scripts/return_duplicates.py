#!/usr/bin/env
"""
Add duplicates back in to a FASTQ file

The output FASTQ file will contain 'cluster_size' copies of each read,
as described by the read-description JSON file
"""

from json import load
from copy import copy
from argparse import ArgumentParser
from itertools import chain
import sys

from Bio.SeqIO import parse, write


def parse_args():
    arg_parser = ArgumentParser()
    arg_parser.add_argument('--input-fastq', help='De-duplicated FASTQ')
    arg_parser.add_argument('--read-description',
                            help="JSON file describing the reads.  Must contain the fields 'id' and 'cluster_size'")
    arg_parser.add_argument('--output-fastq', help='FASTQ file to write duplicates to')
    return arg_parser.parse_args()


def get_id_to_cluster_size(read_description_file):
    """
    Input:
    {
        "rows": [
                    {"id": "@SRR1", "cluster_size": 4},
                    {"id": "@SRR2", "cluster_size": 1}
        ]
    }

    Output:
    {'@SRR1': 4, '@SRR2': 1}

    :param read_description_file: JSON filename
    :return: dict
    """
#    with open(read_description_file) as f:
#        read_dict = load(f)
#        return {row['id']: int(row['cluster_size']) for row in read_dict['rows']}

    read_dict = {}
    with open(read_description_file) as f:
        for line in f:
            line = line.strip()
            if(line.startswith("@") or not line): continue

            _id = line.split('\t')[0]
            size = int(line.split('\t')[1])

            read_dict[_id] = size

    return read_dict


def duplicate(read, i):
    """
    Copy read and label with a new ID based on i

    Eg. i = 7:
    @SRR10 ->
    @SRR10_7

    @SRR3/2 ->
    @SRR3_7/2

    :param read: Bio.SeqIO.SeqRecord
    :param i: int
    :return: Bio.SeqIO.SeqRecord
    """
    new_read = copy(read)

    new_id = '{pre_slash}_{i}'.format(pre_slash=read.id.split('/')[0], i=i) + \
             ('/{}'.format(read.id.split('/')[1]) if '/' in read.id else '')

    # Have to change both the id and description for Bio.SeqIO
    new_read.id = new_id
    new_read.description = new_id

    return new_read


def with_duplicates(read, cluster_size):
    """
    Duplicate read to have cluster_size many reads

    :param read: Bio.SeqIO.SeqRecord
    :param cluster_size: int
    :return: list of duplicate reads
    """
    duplicates = [duplicate(read, i) for i in range(2, cluster_size + 1)]
    return [read] + duplicates


def add_duplicates(reads, id_to_cluster_size):
    """
    Duplicate reads as necessary to have the amount of reads described by id_to_cluster_size

    :param reads: list of reads
    :param id_to_cluster_size: dict
    :return: list of reads
    """
    # Expand each read into a list of duplicates, then merge into 1 list
    return chain.from_iterable([with_duplicates(read, id_to_cluster_size[read.id]) for read in reads])


args = parse_args()

id_to_cluster_size = get_id_to_cluster_size(args.read_description)

de_duplicated_reads = list(parse(args.input_fastq, 'fastq'))
reads = add_duplicates(de_duplicated_reads, id_to_cluster_size)

write(reads, args.output_fastq, 'fastq')
