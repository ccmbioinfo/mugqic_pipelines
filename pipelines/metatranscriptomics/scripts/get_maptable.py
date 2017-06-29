#!/usr/bin/env python

"""
Generates maptable that shows the number of reads used to assemble each contig

Input: {readset}.contigs.json
Output: {readset}.maptable.txt; tab separated
"""

import json
from argparse import ArgumentParser
from Bio import SeqIO

def parse_args():
    arg_parser = ArgumentParser()
    arg_parser.add_argument('--id-file', help='path to pair of mapped reads in txt')
    arg_parser.add_argument('--contig-fasta', help='path to assembled contig in fasta format')
    arg_parser.add_argument('--maptable', help='maptable file to write to in .txt')

    return arg_parser.parse_args()

def get_id_to_num_reads(id_file):
    """
    Reads contigs_id.json file
    and checks how many times same
    contig_id occurs

    Input: {readset}.trinity_id.txt

    Output: dictionary (contig_id as key, # of occurence as value)
    """
    # Dictionary to hold number of reads per contig
    id_to_num_reads = {}
    with open(id_file, 'r+') as f:
        for line in f:
            line = line.strip()
            # Skip a header line
            if(line.startswith("@") or not line):
                continue

            _id = line.split('\t')[1]

            if(_id in id_to_num_reads):
                id_to_num_reads[_id] += 1
                print(_id)
            else:
                id_to_num_reads[_id] = 1

    return id_to_num_reads

def write_maptable(args, id_to_num_reads, id_to_length):
    """
    Write mapping info as txt
    """
    #with open(args.maptable, 'w+') as fh:
    #    json.dump({
    #        'rows': [
    #            {
    #                'id': contig_id, 'num_reads': id_to_num_reads[contig_id],
    #                'length': id_to_length[contig_id]
    #                } for contig_id in id_to_length
    #            ]
    #        }, fh, indent=4)

    with open(args.maptable, 'w+') as fh:
        # Write header
        fh.write("@id\tnum_reads\tlength\n")
        for contig_id in id_to_length:
            fh.write("{_id}\t{num_reads}\t{length}\n".format(_id=contig_id,
                                                            num_reads=id_to_num_reads[contig_id],
                                                            length=id_to_length[contig_id]))

def get_id_to_length(args, id_to_num_reads):
    """
    Reads input fasta file and stores length of
    sequence in a dictionary as value

    Output: dictionary (key=contig_id, value=length of sequence)
    """
    # dictionary to hold sequence length per each contig
    id_to_length = {}
    with open(args.contig_fasta, "rU") as fh:
        for record in SeqIO.parse(fh, "fasta"):
            contig_id = record.id
            id_to_length[contig_id] = len(record.seq)

            if not contig_id in id_to_num_reads:
                id_to_num_reads[contig_id] = 1

    return id_to_length


args = parse_args()

id_to_num_reads = get_id_to_num_reads(args.id_file)
id_to_length = get_id_to_length(args, id_to_num_reads)
write_maptable(args, id_to_num_reads, id_to_length)
