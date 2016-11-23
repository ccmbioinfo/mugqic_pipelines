#!/usr/bin/env python
"""
Output:
<ID><TAB>[2nd field of id-file<TAB>]<sequence length>

id-file is optional, and fills in the 2nd field
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
    def get_id(line):
        return line.strip().split('\t')[0]

    def get_value(line):
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
# for i in (1, 2):
#     id_filename = os.path.join(input_dir, 'cow{}_qual_all_unique_IDs.txt'.format(i))
#     fasta_filename = os.path.join(input_dir, 'cow{}_qual_all_unique.fasta'.format(i))
#
#     out_filename = os.path.join(output_dir, 'cow{}_IDs_length.txt'.format(i))
#
#     reads = dict()
#     with open(id_filename) as id_file:
#         for line in id_file:
#             id, value = line.rstrip().split('\t')
#             reads[id] = value
#
#     print('Number of unique reads: ' + str(len(reads)))
#
#     print('Writing to ' + out_filename)
#     with open(out_filename, 'w+') as out_file:
#         for record in SeqIO.parse(fasta_filename, 'fasta'):
#             out_file.write('\t'.join([record.id, str(reads[record.id]), str(len(record.seq))]) + '\n')
