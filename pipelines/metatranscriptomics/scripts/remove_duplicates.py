#!/usr/bin/env python
"""
Remove duplicates from a fastq file

Record the IDs, cluster size, and the length of each sequence in a *.json file
"""
import json
import re
from argparse import ArgumentParser

from Bio import SeqIO


def parse_args():
    arg_parser = ArgumentParser()
    arg_parser.add_argument('--input-fastq', help='Original fastq before de-duplication')
    arg_parser.add_argument('--unique-fasta', help='Fasta that has been de-duplicated')
    arg_parser.add_argument('--unique-uc', help='UC output from duplicate clustering step')

    arg_parser.add_argument('--output-read-description', help='File to store the IDs of reads that have been removed')
    arg_parser.add_argument('--output-fastq', help='Fastq that has been de-duplicated')
    arg_parser.add_argument('--output-fasta', help='Fasta that has been de-duplicated')

    return arg_parser.parse_args()


def get_unique_ids(uc_file):
    """
    Get the fastq IDs of reads that are unique
    """
    def get_id(line):
        return re.match('^(\S+\t){8}(\S+)', line).group(2)

    return {get_id(line) for line in open(uc_file) if line.startswith('S')}


def get_cluster_id_to_size(fasta):
    """
    ID of a cluster -> # of reads in cluster

    :param fasta: fasta that has been clustered with usearch
    :return: dict
    """
    def get_id(line):
        return re.search('Cluster(\d+);', line).group(1)

    def get_size(line):
        return re.search('size=(\d+)', line).group(1)

    return {get_id(line): get_size(line) for line in open(fasta) if line.startswith('>')}


def get_fastq_id_to_cluster_id(uc_file):
    """
    Fastq ID -> ID of cluster that read is part of

    :param uc_file: *.uc from usearch
    :return: dict
    """
    def fastq_id(line):
        return re.match('^(\S+\t){8}(\S+)', line).group(2)

    def cluster_id(line):
        return re.match('^(\S+\t){2}(\S+)', line).group(2)

    return {fastq_id(line): cluster_id(line) for line in open(uc_file) if line.startswith('S')}


def get_unique_reads(fastq, uc_file):
    """
    Subset the reads from fastq to remove duplicates

    :param fastq: original fastq with duplicates
    :param uc_file: *.uc file
    :return: list
    """
    unique_ids = get_unique_ids(uc_file)
    return [read for read in SeqIO.parse(fastq, 'fastq') if read.id in unique_ids]


def assign_cluster_size(unique_reads, uc_file, unique_fasta):
    """
    Set the 'cluster_size' attribute on all unique reads

    :param unique_reads: iterable of reads
    :param uc_file: *.uc file
    :param unique_fasta: output from usearch that contains the number of duplicates for each read
    """
    fastq_id_to_cluster_id = get_fastq_id_to_cluster_id(uc_file)
    cluster_id_to_size = get_cluster_id_to_size(unique_fasta)

    for read in unique_reads:
        read.cluster_size = cluster_id_to_size[fastq_id_to_cluster_id[read.id]]


def write_id_to_cluster_size(reads, id_file):
    """
    Write a table in json containing fastq ID and cluster size
    Eg:
    {"rows": [
        {"id": "SRR1", "cluster_size": 3}
        {"id": "SRR23", "cluster_size": 1}
    ]}
    """
    id_to_cluster_size = {'rows':
                              [{'id': read.id, 'cluster_size': read.cluster_size, 'length': len(read.seq)} for read in reads]
                          }

    with open(id_file, 'w+') as f:
        json.dump(id_to_cluster_size, f, indent=4)


def write_unique_reads(reads, file, format):
    """
    Write reads to file according to formatt

    :param reads: iterable of reads
    :param file: *.fastq file or *.fasta
    :param format: 'fastq' or 'fasta'
    :return:
    """
    with open(file, 'w+') as out:
        SeqIO.write(reads, out, format)


args = parse_args()

# Get the unique reads and determine how many there are in each cluster
unique_reads = get_unique_reads(args.input_fastq, args.unique_uc)
assign_cluster_size(unique_reads, args.unique_uc, args.unique_fasta)

# Keep track of how many duplicates there are of each read
write_id_to_cluster_size(unique_reads, args.output_read_description)

# Write out a *.fasta and *.fastq file containing only the unique reads
write_unique_reads(unique_reads, args.output_fastq, 'fastq')
write_unique_reads(unique_reads, args.output_fasta, 'fasta')
