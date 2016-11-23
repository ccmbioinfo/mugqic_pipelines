#!/usr/bin/env python
"""
Split fasta or fastq into 2 output files, based on the given IDs in id_filename
See --help
"""
import argparse
from Bio import SeqIO


# Parse args
arg_parser = argparse.ArgumentParser()
arg_parser.add_argument('--fastq')
arg_parser.add_argument('--fasta')
arg_parser.add_argument('--id-file')
arg_parser.add_argument('--included', help='Output filename for reads that are in the IDs')
arg_parser.add_argument('--excluded', help='Output filename for reads that are not in the IDs')
args = arg_parser.parse_args()


# Get the IDs from the id_file
ids = {line.split('\t')[0] for line in open(args.id_file)}

input, format = (args.fastq, 'fastq') if args.fastq else (args.fasta, 'fasta')
all_reads = set(SeqIO.parse(input, format))

included_reads = {read for read in all_reads if read.id in ids}
excluded_reads = all_reads - included_reads

SeqIO.write(included_reads, args.included, format)
SeqIO.write(excluded_reads, args.excluded, format)
# Write to either included_filename or excluded filename depending on whether id is in ids
# with open(args.included_filename, 'w+') as included_file, open(args.excluded_filename, 'w+') as excluded_file:
#     for read in SeqIO.parse(args.fastq_filename, 'fastq'):
#         included_file.write(read.format('fastq')) if read.id in ids else excluded_file.write(read.format('fastq'))
