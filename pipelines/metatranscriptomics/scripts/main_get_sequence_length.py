#!/usr/bin/env python
"""
Read a fastq or fasta and file and write out an ID file containing:
<ID><TAB>[2nd field of id-file<TAB>]<sequence length>

Additional information may be specified by an id-file, which fills in the second field above
"""
import argparse
import os
import json
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
  #  if args.id_file:
  #      def write_read(read):
  #           return '{id}\t{value}\t{len}\n'.format(id=read.id, value=read.value, len=len(read.seq))
  #  else:
  #      def write_read(read):
  #           return '{id}\t{len}\n'.format(id=read.id, len=len(read.seq))

    def write_read(read):
        with_id = '{id}\t{value}\t{len}\n'.format(id=read.id, value=read.value, len=len(read.seq))
        without_id = '{id}\t{len}\n'.format(id=read.id, len=len(read.seq))

        return with_id if args.id_file else without_id

    if('nr_all_sub' in output_file):
        with open(output_file, 'w+') as out:
            # Write header line
            out.write("@id\tdescription\tlength\n")
            for read in reads:
                out.write("{_id}\t{description}\t{length}\n".format(
                                                                    _id=read.id,
                                                                    description=read.description,
                                                                    length=len(read.seq)))
    elif('microbial_cds' in output_file):
        with open(output_file, 'w+') as out:
            # Write header line
            out.write("@id\tlength\n")
            for read in reads:
                out.write("{_id}\t{length}\n".format(_id=read.id,
                                                        length=len(read.seq)))

    else:
        with open(output_file, 'w+') as out:
            # Write header line
            out.write("@id\tlength\n")
            for read in reads:
                out.write("{_id}\t{length}\n".format(_id=read.id,
                                                        length=len(read.seq)))


args = parse_args()
input, format = (args.fastq, 'fastq') if args.fastq else (args.fasta, 'fasta')

reads = list(SeqIO.parse(input, format))

if args.id_file:
    reads = parse_id_file(args.id_file, reads)

write_output(reads, args.output)
