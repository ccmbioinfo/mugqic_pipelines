"""
Reads accession2taxID files and generate microbial_taxID_all.txt database

Input:
    accession2tax files downloaded from NCBI Taxonomy database
    ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/

    ID_file generated from extract_gi_species_from_fasta.py

Output:
    Tab delimited txt file with gi|gi#|ref|accID|, gi#, Species Name, taxID
"""

from argparse import ArgumentParser

def parse_args():
    arg_parser = ArgumentParser()
    arg_parser.add_argument('--tax-gi', nargs='+', help='accession2tax files')
    arg_parser.add_argument('--id-file', help='File containing gi|#|ref|Acc|, gi, species name')
    arg_parser.add_argument('--out', help='output name')

    return arg_parser.parse_args()

def read_tax_file(args):
    tax_dict = {}

    for tax_gi_file in args.tax_gi:
        with open(tax_gi_file) as fh:
            # Skip header
            next(fh)
            for line in fh:
                line = line.strip()
                if(not line): continue
                _gi = int(line.split('\t')[3])
                _tax = line.split('\t')[2]

                tax_dict[_gi] = _tax

        print("{fh}: Processing Done!".format(fh=tax_gi_file))

    return tax_dict

def append_tax_to_id(args, tax_dict):
    result = []

    with open(args.id_file) as fh:
        for line in fh:
            line = line.strip()
            if(not line): continue

            gi = int(line.split('\t')[1])

            if(gi in tax_dict):
                result.append('\t'.join((line, tax_dict[gi])))

    return result

def write(args, result):
    with open(args.out, 'w+') as fh:
        for i in range(0,len(result)):
            fh.write("{info}\n".format(info=result[i]))

args = parse_args()

tax_dict = read_tax_file(args)

result = append_tax_to_id(args, tax_dict)

write(args,result)

