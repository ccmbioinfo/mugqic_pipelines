"""
Merge multiple ID files into one and remove duplicate IDs.
"""

from argparse import ArgumentParser
import glob
import os

def parse_args():
    arg_parser = ArgumentParser()
    arg_parser.add_argument('--alignment', nargs='*', type=str, help='user '
            'specified pattern to merge certain ID files')
    arg_parser.add_argument('--location', help='path to directory with ID files')
    arg_parser.add_argument('--unique', help='Merged ID file with no duplicates')

    return arg_parser.parse_args()

def merge(args):
    """
    Merge user specified ID files and store it as a list
    """
    result = []

    file_to_parse = []

    location = os.path.abspath(args.location)
    for pattern in args.alignment:
        file_to_parse.extend(glob.glob(''.join((location,'/','*',pattern,'*','.txt'))))

    for f in file_to_parse:
        with open(f, 'r') as fh:
            for line in fh:
                line = line.strip()
                # Skip if header or empty line
                if(line.startswith("@") or not line): continue

                result.append(line)


    #with open(args.merged, 'w+') as out:
    #    json.dump({
    #        'rows': [
    #            {key: data[key] for key in data.keys() } for data in result
    #            ]
    #        }, out, indent=4)

    return result

def write_unique_hit(merged):
    def get_hit(line):
        return line.split('\t')[1]

    unique = set( get_hit(line) for line in merged )

    with open(args.unique,'w+') as out:
        for hit in unique:
            out.write(hit + '\n')

args = parse_args()

merged = merge(args)
write_unique_hit(merged)
