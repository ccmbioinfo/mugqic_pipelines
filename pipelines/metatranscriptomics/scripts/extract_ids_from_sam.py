#!/usr/bin/env
"""
Read a sam file and extract the IDs to a JSON file
"""

import json
from argparse import ArgumentParser


HEADER_SYMBOL = '@'


def parse_args():
    arg_parser = ArgumentParser(description=__doc__)
    arg_parser.add_argument('--sam', help='Output from samtools')
    arg_parser.add_argument('--id-file', help='JSON file containing read IDs')
    return arg_parser.parse_args()


def get_ids(sam):
    """
    Get the IDs from the given sam file

    Example input:
    @SQ                ...
    @SQ                ...
    @PG                ...
    SRR5       77      ...
    SRR6       141     ...
    SRR7       77      ...

    Example output:
    {'SRR5, 'SRR6', 'SRR7'}

    :param sam: *.sam filename
    :return: set
    """
    def get_id(line):
        return line.split('\t')[0]

    return {get_id(line) for line in open(sam) if not line.startswith(HEADER_SYMBOL)}


def write_ids(ids, file):
    """
    Write ids to a JSON file

    Example input:
    {'@SRR1', '@SRR2'}

    Example JSON output:
    {
        "rows": [
            {"id": "@SRR1"},
            {"id": "@SRR2"}
        ]
    }

    :param ids: iterable of IDs
    :param file: JSON filename
    """
    with open(file, 'w+') as f:
        json.dump({
            'rows': [
                {'id': id} for id in ids
            ]
        }, f, indent=4)


args = parse_args()
ids = get_ids(args.sam)
write_ids(ids, args.id_file)
