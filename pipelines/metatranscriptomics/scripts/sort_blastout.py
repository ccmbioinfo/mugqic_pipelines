"""
Sort alignment output files
"""

from os import system
from argparse import ArgumentParser
import subprocess

def parse_args():
    arg_parser = ArgumentParser()
    arg_parser.add_argument('--infile', help='output from alignment tools (bwa, blat, etc.)')
    arg_parser.add_argument('--sort', help='file to write the result to')
    arg_parser.add_argument('--cutoff', help='cutoff threshold')

    return arg_parser.parse_args()

def write_file(args, result):
    with open(args.sort, 'w+') as fh:
        for val in result:
            fh.write(str(val) + '\n')

def read_file(args):
    id_dict = {}
    result = []
    with open(args.infile, 'r') as fh:
        for line in fh:
            head = line.split('\t')[0].split('/')

            if(head[0] not in id_dict):
                id_dict[head[0]] = 1

                hit_pairs = subprocess.check_output(['grep',head[0],args.infile]) \
                                    .strip() \
                                    .split('\n')

                length = len(hit_pairs)

                if(length == 1):
                    result.append(hit_pairs[0])
                else:
                    pairs = {}
                    scores = {}

                    for k in range(length):
                        pair = hit_pairs[k].split('\t')
                        pairs[k] = hit_pairs[k]
                        scores[k] = pair[11]

                    sorted_scores_keys = sorted(scores, reverse=True, key=scores.get)
                    for i in range(0, min(length, int(args.cutoff))):
                        result.append(pairs[sorted_scores_keys[i]])
    return result

args = parse_args()
result = read_file(args)
write_file(args, result)
