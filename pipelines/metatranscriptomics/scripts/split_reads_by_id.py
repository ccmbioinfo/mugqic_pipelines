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
print('Total number of reads: {}'.format(len(all_reads)))

included_reads = {read for read in all_reads if read.id in ids}
excluded_reads = all_reads - included_reads

print('Number of included reads: {}'.format(len(included_reads)))
print('Number of excluded reads: {}'.format(len(excluded_reads)))

SeqIO.write(included_reads, args.included, format)
SeqIO.write(excluded_reads, args.excluded, format)
