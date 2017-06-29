#!/usr/bin/env python


import json
from argparse import ArgumentParser

DEFAULT_READ_LENGTH = 100


def parse_args():
    arg_parser = ArgumentParser()
    arg_parser.add_argument('--read-description',
                            help="TXT file describing the unique reads.  Contains 'id' and 'length'")
    arg_parser.add_argument('--infernalout', help='Output from cmscan')
    arg_parser.add_argument('--apply-cutoff', help='Whether to apply cutoff values', action='store_true')
    arg_parser.add_argument('--max-evalue', help='The maximum e-value to accept a read as being rRNA')
    arg_parser.add_argument('--min-percent-identity', help='Minimal percentage identity to accept a read as being rRNA')

    arg_parser.add_argument('--out-ids', help='Path to write IDs that are rRNA in JSON')
    return arg_parser.parse_args()


def get_id_to_length(read_description_file):
    """
    Return a dict mapping from 'id' to 'length'

    Example JSON file:
    {
        "rows": [{"id": "@SRR1", "length": 100},
                 {"id": "@SRR2", "length": 98}]
    }

    Example output:
    {'@SRR1': 100, '@SRR2': 98}

    :param read_description_file: JSON file
    :return: dict
    """
#    with open(read_description_file) as f:
#        read_dict = json.load(f)
#        return {row['id']: int(row['length']) for row in read_dict['rows']}
    read_dict = {}
    with open(read_description_file) as f:
        for line in f:
            line = line.strip()
            if(line.startswith("@") or not line): continue

            _id = line.split('\t')[0]
            length = int(line.split('\t')[2])

            read_dict[_id] = length

    return read_dict


class InfernalParser:
    """
    Bunch of staticmethod's for parsing the output of infernal

    Example:
    InfernalParser.get_id("my line")
    InfernalParser.get_evalue("my other line")
    """
    @staticmethod
    def get_id(line):
        return line.split()[2]

    @staticmethod
    def get_seq_from(line):
        return int(line.split()[7])

    @staticmethod
    def get_seq_to(line):
        return int(line.split()[8])

    @staticmethod
    def get_evalue(line):
        return float(line.split()[15])


def get_percent_identity(seq_from, seq_to, seq_length):
    return 100 * abs(seq_to - seq_from + 1) / seq_length


def passes_cutoff(evalue, percent_identity, args):
    return evalue <= float(args.max_evalue) and percent_identity >= float(args.min_percent_identity)


def parse_infernalout(infernalout, id_to_length, args):
    """
    Read the infernalout file and determine which read IDs are rRNA

    :param infernalout: table file from cmscan
    :param id_to_length: dict
    :return: set of str
    """
    rrna_ids = set()

    with open(infernalout) as f:
        for line in f:
            # Skip comment lines
            if line.startswith('#'): continue

            id = InfernalParser.get_id(line)
            seq_from = InfernalParser.get_seq_from(line)
            seq_to = InfernalParser.get_seq_to(line)
            evalue = InfernalParser.get_evalue(line)

            if id not in rrna_ids:
                if id not in id_to_length:
                    id_to_length[id] = DEFAULT_READ_LENGTH

                percent_identity = get_percent_identity(seq_from, seq_to, id_to_length[id])

                if not args.apply_cutoff or passes_cutoff(evalue, percent_identity, args):
                    rrna_ids.add(id)

    return rrna_ids


def write_rrna_ids(ids, output_ids):
    """
    Write IDs to JSON file

    Example input:
    {'@SRR1', '@SRR5'}

    Example JSON output:
    {
        "rows": [{"id": "@SRR1"},
                 {"id": "@SRR5"}]
    }

    :param ids: iterable of str
    :param output_ids: JSON filename
    """
#    with open(output_ids, 'w+') as f:
#        json.dump({
#            'rows': [{'id': id} for id in ids]
#        }, f, indent=4)

    with open(output_ids, 'w+') as f:
        f.write("@id\n")
        for _id in ids:
            f.write("{_id}\n".format(_id=_id))


args = parse_args()

id_to_length = get_id_to_length(args.read_description)
rrna_ids = parse_infernalout(args.infernalout, id_to_length, args)
write_rrna_ids(rrna_ids, args.out_ids)
