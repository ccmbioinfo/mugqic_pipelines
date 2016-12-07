#!/usr/bin/env python


import json
from argparse import ArgumentParser

DEFAULT_READ_LENGTH = 100


def parse_args():
    arg_parser = ArgumentParser()
    arg_parser.add_argument('--read-description',
                            help="JSON file describing the unique reads.  Contains 'id', 'cluster_size', and 'length'")
    arg_parser.add_argument('--infernalout', help='Output from cmscan')
    arg_parser.add_argument('--apply-cutoff', help='Whether to apply cutoff values', action='store_true')
    arg_parser.add_argument('--max-evalue', help='The maximum e-value to accept a read as being rRNA')
    arg_parser.add_argument('--min-percent-identity', help='Minimal percentage identity to accept a read as being rRNA')

    arg_parser.add_argument('--out-ids', help='Path to write IDs that are rRNA in JSON')
    return arg_parser.parse_args()


def parse_read_description(read_description_file):
    """
    Index a dict mapping from 'id' -> {'cluster_size': <cluster_size>, 'length': <length>}

    Example JSON file:
    {
        "rows": [{"id": "@SRR1", "cluster_size": 2, "length": 100},
                 {"id": "@SRR2", "cluster_size": 1, "length": 98}]
    }

    Example output:
    {'@SRR1': {'cluster_size': 2, 'length': 100},
     '@SRR2': {'cluster_size': 1, 'length': 98}}

    :param read_description_file: JSON file
    :return: dict
    """
    with open(read_description_file) as f:
        read_dict = json.load(f)
        return {row['id']: {'length': row['length'], 'cluster_size': row['cluster_size']} for row in read_dict['rows']}


class InfernalFileParser:
    @staticmethod
    def get_id(line):
        return line.split()[1]

    @staticmethod
    def get_seq_from(line):
        return line.split()[7]

    @staticmethod
    def get_seq_to(line):
        return line.split()[8]

    @staticmethod
    def get_score(line):
        return line.split()[14]

    @staticmethod
    def get_evalue(line):
        return line.split()[15]


def get_percent_identity(seq_from, seq_to, seq_length):
    # TODO: round to 2 decimal places
    return 100 * abs(seq_to - seq_from + 1) / seq_length


def passes_cutoff(evalue, percent_identity, args):
    return evalue <= float(args.max_evalue) and percent_identity >= float(args.min_percent_identity)


def parse_infernalout(infernalout, read_description, args):
    """
    Read the infernalout file and determine which read IDs are rRNA

    :param infernalout: table file from cmscan
    :param read_description: dict
    :return: set of str
    """
    rrna_ids = set()

    with open(infernalout) as f:
        for line in f:
            id = InfernalFileParser.get_id(line)
            seq_from = InfernalFileParser.get_seq_from(line)
            seq_to = InfernalFileParser.get_seq_to(line)
            score = InfernalFileParser.get_score(line)
            evalue = InfernalFileParser.get_evalue(line)

            if id not in rrna_ids:
                if not hasattr(read_description[id], 'length'):
                    read_description[id]['length'] = DEFAULT_READ_LENGTH

                percent_identity = get_percent_identity(seq_from, seq_to, int(read_description[id]['length']))

                if not args.apply_cutoff or passes_cutoff(percent_identity, evalue, args):
                    rrna_ids.add(id)


def write_rrna_ids(rrna_ids, output_ids):
    with open(output_ids, 'w+') as f:
        json.dump({
            'rows': [{'id': id} for id in rrna_ids]
        }, f)

args = parse_args()

read_description = parse_read_description(args.read_description)
rrna_ids = parse_infernalout(args.infernalout, read_description, args)
write_rrna_ids(rrna_ids, args.out_ids)
