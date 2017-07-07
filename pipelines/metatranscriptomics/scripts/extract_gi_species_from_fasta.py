"""
Extract geneID, giID, speicies_name, and taxonID from microbial database

Input:
    Microbial database downloaded from NCBI
    ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_refseq/Bacteria/all.ffn.tar.gz

Output:
    Tab delimited txt file with gi|gi#|ref|accID|, gi#, Species Name

Example Command:
    python extract_gi_species_from_fasta.py --fasta microbial_all.fasta --out id_file
"""

from argparse import ArgumentParser

def parse_args():
    arg_parser = ArgumentParser()
    arg_parser.add_argument('--fasta', help='Fasta file to extract information from')
    arg_parser.add_argument('--out', help='output name')

    return arg_parser.parse_args()

def extract_data(args):
    data = []
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
