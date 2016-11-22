#!/usr/bin/env python
import re
import os
from Bio import SeqIO


def get_cluster_size(record):
    # print(record.description)
    return re.search('size=(\d+)', record.description).group(1)


def get_unique_ids(uc_file):
    """
    Get the fastq IDs of reads that are unique
    """
    # ids = set()
    # with open(file) as f:
    #     for line in f:
    #         if line.startswith('S'):
    #             ids.add(re.match('^(\S+\t){8}(\S+)', line).group(2))
    # return ids
    def get_id(line):
        return re.match('^(\S+\t){8}(\S+)', line).group(2)

    return {get_id(line) for line in open(uc_file) if line.startswith('S')}


def get_cluster_id_to_size(fasta):
    cluster_id_to_size = dict()

    with open(fasta) as fasta:
        for line in fasta:
            match = re.search('Cluster(\d+);size=(\d+)', line)
            if match is not None:
                cluster_id_to_size[match.group(1)] = match.group(2)

    return cluster_id_to_size


def get_fastq_id_to_cluster_id(uc_file):
    def fastq_id(line):
        return re.match('^(\S+\t){8}(\S+)', line).group(2)

    def cluster_id(line):
        return re.match('^(\S+\t){2}(\S+)', line).group(2)

    return {fastq_id(line): cluster_id(line) for line in open(uc_file) if line.startswith('S')}


for i in (1, 2):
    original_fastq = 'flash/cow{}_qual_all.fastq'.format(i)
    unique_fasta = 'cluster_duplicates/cow{}_qual_all_unique.fasta'.format(i)
    unique_uc = 'cluster_duplicates/cow{}_qual_all_unique.uc'.format(i)

    output_id_file = 'remove_duplicates/cow{}_qual_all_unique_IDs.txt'.format(i)
    output_fastq = 'remove_duplicates/cow{}_qual_all_unique.fastq'.format(i)
    output_fasta = 'remove_duplicates/cow{}_qual_all_unique.fasta'.format(i)

    unique_ids = get_unique_ids(unique_uc)
    # unique_ids = {recrd.id for recrd in SeqIO.parse('remove_duplicates/cow{}_qual_all_unique.fasta', 'fasta')}
    # Construct a dictionary indexed by sequence
    # unique_reads = {read.seq: read for read in SeqIO.parse(original_fastq, 'fastq') if read.id in unique_ids}
    unique_reads = [read for read in SeqIO.parse(original_fastq, 'fastq') if read.id in unique_ids]
    print('Number of unique reads for cow{}: {}'.format(i, len(unique_reads)))

    fastq_id_to_cluster_id = get_fastq_id_to_cluster_id(unique_uc)
    cluster_id_to_size = get_cluster_id_to_size(unique_fasta)

    for read in unique_reads:
        read.cluster_size = cluster_id_to_size[fastq_id_to_cluster_id[read.id]]

    # Write the *_qual_all_unique_IDs.txt file as <ID><TAB><number of duplicates>
    with open(output_id_file, 'w+') as f:
        for read in unique_reads:
            f.write('{id}\t{cluster_size}\n'.format(id=read.id, cluster_size=read.cluster_size))

    # Overwrite the unique fastq file
    with open(output_fastq, 'w+') as out:
        SeqIO.write(unique_reads, out, 'fastq')

    # Overwrite the unique fasta file
    with open(output_fasta, 'w+') as out:
        SeqIO.write(unique_reads, out, 'fasta')
