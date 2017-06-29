"""
Extract geneID, giID, speicies_name, and taxonID from microbial database
"""

from argparse import ArgumentParser
import re

def parse_args():
    arg_parser = ArgumentParser()
    arg_parser.add_argument('--fasta', help='Fasta file to extract information from')
    arg_parser.add_argument('--out', help='output format')

    return arg_parser.parse_args()

def extract_data(args):
    data = []
    prev = ''
    with open(args.fasta, 'r') as fh:
        for line in fh:
            line = line.strip()
            giID = line.split('|')[1]

            fullID = line.split(' ')[0]
            fullID_formatted = fullID.replace('>', '')
            species = line.replace(fullID, '').strip()
#            species = re.search(':.*?([A-Z].*?)($|chromosome|,)',line).group(1) \
#                        .strip()

            data.append('\t'.join((fullID_formatted, giID, species)))

    return data

def write_file(args, data):
    with open(args.out, 'w+') as fh:
        for i in range(0, len(data)):
            fh.write('{data}\n'.format(data=data[i]))

args = parse_args()

data = extract_data(args)

write_file(args, data)
