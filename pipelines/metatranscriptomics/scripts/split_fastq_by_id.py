#!/usr/bin/env python
"""
Split fastq_filename into 2 output files, based on the given IDs in id_filename
See --help
"""

import argparse
from Bio import SeqIO


arg_parser = argparse.ArgumentParser()
arg_parser.add_argument('fastq_filename')
arg_parser.add_argument('id_filename')
arg_parser.add_argument('included_filename')
arg_parser.add_argument('excluded_filename')
args = arg_parser.parse_args()


# Get the IDs from the id_file
with open(args.id_filename) as id_file:
    ids = {line.split('\t')[0] for line in id_file}

# print('Total number of IDs: {}'.format(len(ids)))


# Write to either included_filename or excluded filename depending on whether id is in ids
with open(args.included_filename, 'w+') as included_file, open(args.excluded_filename, 'w+') as excluded_file:
    for read in SeqIO.parse(args.fastq_filename, 'fastq'):
        included_file.write(read.format('fastq')) if read.id in ids else excluded_file.write(read.format('fastq'))
