#!/usr/bin/env python
"""
Read a fastq or fasta and file and write out an ID file containing:
<ID><TAB>[2nd field of id-file<TAB>]<sequence length>

Additional information may be specified by an id-file, which fills in the second field above
"""
import argparse
import os
from Bio import SeqIO


def parse_args():
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument('--fastq')
    arg_parser.add_argument('--fasta')
    arg_parser.add_argument('--id-file')
    arg_parser.add_argument('--output')
    return arg_parser.parse_args()


def parse_id_file(id_file, reads):
    """
    Write the "value" field of read based on id_file
    :param reads: list of reads from SeqIO
    :return:
    """
    def get_id(line):
        """
        Get the ID from a line
        """
        return line.strip().split('\t')[0]

    def get_value(line):
        """
        Get the 2nd value from a line
        """
        return line.strip().split('\t')[1]

    id_to_value = {get_id(line): get_value(line) for line in open(id_file)}

    # Set the "value" field of each read
    for read in reads:
        read.value = id_to_value[read.id]

    return reads


def write_output(reads, output_file):
    if args.id_file:
        def write_read(read):
            return '{id}\t{value}\t{len}\n'.format(id=read.id, value=read.value, len=len(read.seq))
    else:
        def write_read(read):
            return '{id}\t{len}\n'.format(id=read.id, len=len(read.seq))

    with open(output_file, 'w+') as out:
        for read in reads:
            out.write(write_read(read))


args = parse_args()
input, format = (args.fastq, 'fastq') if args.fastq else (args.fasta, 'fasta')

reads = list(SeqIO.parse(input, format))

if args.id_file:
    reads = parse_id_file(args.id_file, reads)

write_output(reads, args.output)
