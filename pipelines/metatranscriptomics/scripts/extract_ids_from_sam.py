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

def get_pairs(sam):
    """
    If the input is host.bwaout, get the IDs from the given sam file

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

    If the input is contigs.bwaout, get putative mRNA IDs and its contigs ID

    Example contigs.bwaout input
    SRR594215.107462    145    TRINITY_DN96_c0_g1_i1    ...
    SRR594215.108123    163    TRINITY_DN11_c0_g1_i1    ...
    SRR594215.118549    177    TRINITY_DN262_c0_g1_i1    ...
    SRR594215.123113    73    TRINITY_DN49_c0_g2_i1    ...
    SRR594215.126943    137    TRINITY_DN50_c1_g3_i1    ...

    Example output:
    {SRR594215.107462: TRINITY_DN96_c0_g1_i1}

    return dictionary
    """
    def get_id(line):
        return line.split('\t')[0]
    def get_hit(line):
        return line.split('\t')[2]

    return {get_id(line) : get_hit(line) for line in open(sam) if not
            line.startswith(HEADER_SYMBOL)}

def write_id_file(pairs, file):
    """
    If writing host_id.json

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

    If writing contigs_id.json

    Example input
    {SRR594215.107462: TRINITY_DN96_c0_g1_i1}

    Example output
    {
        "rows": [
            {
                "contig_id": "TRINITY_DN43_c0_g7_i1",
                "mRNA_id": "SRR594215.3009418_2"
            },
        ]
    }
    """
#    if "contigs.json" in file:
#        with open(file, 'w+') as f:
#            json.dump({
#                'rows': [
#                    {'mRNA_id': id, 'contig_id': pairs[id]} for id in pairs
#                ]
#            }, f, indent=4)
#    else:
#        with open(file, 'w+') as f:
#            json.dump({
#                'rows': [
#                    {'id': id, 'hit': pairs[id]} for id in pairs
#                ]
#            }, f, indent=4)

    with open(file, 'w+') as f:
        f.write("@id\thit\n")
        for _id in pairs:
            f.write("{_id}\t{hit}\n".format(_id=_id,
                                            hit=pairs[_id]))

args = parse_args()
pairs = get_pairs(args.sam)
write_id_file(pairs, args.id_file)
