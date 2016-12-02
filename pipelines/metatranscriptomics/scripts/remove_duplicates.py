#!/usr/bin/env python
"""
Input:
Original fastq          - *.{1,2}.qual_all.fastq
Only unique reads       - *.{1,2}.qual_all_unique.fasta
UC file of unique reads - *.{1,2}.qual_all_unique.uc

Output:
ID to # of duplicates                    - *.{1,2}.qual_all_unique_IDs.txt
De-duplicated fastq                      - *.{1,2}.qual_all_unique.fastq
De-duplicated fasta (with corrected IDs) - *.{1,2}.qual_all_unique.fasta
"""
import re
from argparse import ArgumentParser

from Bio import SeqIO


def parse_args():
    arg_parser = ArgumentParser()
    arg_parser.add_argument('--input-fastq', help='Original fastq before de-duplication')
    arg_parser.add_argument('--unique-fasta', help='Fasta that has been de-duplicated')
    arg_parser.add_argument('--unique-uc', help='UC output from duplicate clustering step')

    arg_parser.add_argument('--output-ids', help='File to store the IDs of reads that have been removed')
    arg_parser.add_argument('--output-fastq', help='Fastq that has been de-duplicated')
    arg_parser.add_argument('--output-fasta', help='Fasta that has been de-duplicated')

    return arg_parser.parse_args()


def get_cluster_size(record):
    # print(record.description)
    return re.search('size=(\d+)', record.description).group(1)


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

    return {get_id(line): get_size(line) for line in open(fasta)}


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
    unique_ids = get_unique_ids(uc_file)
    return [read for read in SeqIO.parse(fastq, 'fastq') if read.id in unique_ids]


def assign_cluster_size(reads, uc_file, unique_fasta):
    fastq_id_to_cluster_id = get_fastq_id_to_cluster_id(uc_file)
    cluster_id_to_size = get_cluster_id_to_size(unique_fasta)

    for read in reads:
        read.cluster_size = cluster_id_to_size[fastq_id_to_cluster_id[read.id]]


def write_ids(reads, id_file):
    # Write the ID file as <ID><TAB><number of duplicates>
    with open(id_file, 'w+') as f:
        for read in reads:
            f.write('{id}\t{cluster_size}\n'.format(id=read.id, cluster_size=read.cluster_size))


def write_unique_reads(reads, file, format):
    with open(file, 'w+') as out:
        SeqIO.write(reads, out, format)


args = parse_args()

unique_reads = get_unique_reads(args.original_fastq, args.unique_uc)
assign_cluster_size(unique_reads, args.unique_uc, args.unique_fasta)

write_ids(unique_reads, args.output_ids)
write_unique_reads(unique_reads, args.unique_fastq, 'fastq')
write_unique_reads(unique_reads, args.unique_fasta, 'fasta')
