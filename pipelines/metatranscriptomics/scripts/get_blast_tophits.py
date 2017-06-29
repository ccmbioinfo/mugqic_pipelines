"""
Extract top hits (Similar to get_blast_1tophit.py but there's slight difference)
"""

import json
from argparse import ArgumentParser

def parse_args():
    arg_parser = ArgumentParser()
    arg_parser.add_argument('--length', help='length file')
    arg_parser.add_argument('--infile', help='input file')
    arg_parser.add_argument('--outfile', help='file to output data')
#    arg_parser.add_argument('--json', help='json file to write id for future processing')
    arg_parser.add_argument('--id-file', help='txt file to write id for future processing')
    arg_parser.add_argument('--cutoff-type')
    arg_parser.add_argument('--cutoff0')
    arg_parser.add_argument('--cutoff1')
    arg_parser.add_argument('--cutoff2')
    arg_parser.add_argument('--cutoff3')
    arg_parser.add_argument('--diamond', type=bool, default=False)

    return arg_parser.parse_args()

def check_cutoff_type(args):
    if(args.cutoff_type == 0):
        args.cutoff0 = 0
        args.cutoff1 = 0
        args.cutoff2 = 0
        args.cutoff3 = 0

def read_length_file(length_file):
    reads = {}
    with open(length_file) as fh:
        for line in fh:
            line = line.strip()
            # Skip if empty line or header
            if(line.startswith("@") or not line): continue

            _id = line.split('\t')[0]
            if("mappingtable" in length_file):
                length = int(line.split('\t')[2])
            else:
                length = int(line.split('\t')[1])

            reads[_id] = length

    return reads

def read_infile(args, reads):
    pairs = {}
    hits = {}
    line1 = 'na'
    max_score = 0
    k = 0

    with open(args.infile, 'r') as fh:
        for my_line in fh:
            my_line = my_line.strip()
            line = my_line.split('\t')

            dbID = line[0]

            if('nr_sub' in args.infile and dbID not in reads):
                dbID = ''.join(line[0],',')
            if('nr_sub' in args.infile and dbID not in reads):
                dbID = ''.join(line[0],';')

            if(line[0] == line1):
                k += 1
            else:
                k = 0

            if(k == 0):
                line1 = dbID
                max_score = float(line[11])

                hitgiID = line[1]

                pident = float(line[2])
                pident = '{:.2f}'.format(pident)

                evalue = line[10]

                if(dbID not in reads):
                    reads[dbID] = 100
                if(args.diamond):
                    poverlap = 100 * (3 * int(line[3])) / reads[dbID]
                else:
                    poverlap = 100 * int(line[3]) / reads[dbID]

                poverlap = '{:.2f}'.format(float(poverlap))

                if(reads[dbID] >= int(args.cutoff0)):
                    if(max_score >= int(args.cutoff3)):
                        pairs[dbID] = '\t'.join(str(x) for x in (dbID, hitgiID, pident, poverlap, evalue, max_score))
                        hits[dbID] = hitgiID
                    else:
                        max_score += 1
                else:
                    if((pident >= int(args.cutoff1)) and
                    (poverlap >= int(args.cutoff2))):
                        pairs[dbID] = '\t'.join(str(x) for x in (dbID, hitgiID, pident, poverlap, evalue, max_score))
                        hits[dbID] = hitgiID
                    else:
                        max_score += 1
            else:
                if(line[11] == max_score):
                    hitgiID = line[1]

                    pident = float(line[2])
                    pident = '{:.2f}'.format(pident)

                    evalue = line[10]

                    if('diamond' in args.infile):
                        poverlap = 100 * (3 * int(line[3])) / reads[dbID]
                    else:
                        poverlap = 100 * int(line[3]) / reads[dbID]

                    poverlap = '{:.2f}'.format(float(poverlap))

                    if(reads[dbID] >= int(args.cutoff0)):
                        if(max_score >= int(args.cutoff3)):
                            pairs[dbID] = '\t'.join(str(x) for x in (dbID, hitgiID, pident, poverlap, evalue, max_score))
                            hits[dbID] = hitgiID
                    else:
                        if((pident >= int(args.cutoff1)) and
                        (poverlap >= int(args.cutoff2))):
                            pairs[dbID] = '\t'.join(str(x) for x in (dbID, hitgiID, pident, poverlap, evalue, max_score))
                            hits[dbID] = hitgiID

    return pairs, hits

def write_file(args, pairs):
    with open(args.outfile, 'w+') as fh:
        for key in pairs:
            fh.write(pairs[key] + '\n')

def write_id_file(args, hits):
    with open(args.id_file, 'w+') as fh:
        fh.write("@id\thit\n")
        for dbID in hits:
            fh.write('{_id}\t{_hit}\n'.format(_id=dbID,
                                                _hit=hits[dbID]))

# Deprecated
#def write_json(args, hits):
#    with open(args.json, 'w+') as fh:
#        json.dump({
#            'rows': [
#                {
#                    'id': dbID, 'hit': hits[dbID]
#                    } for dbID in hits
#                ]
#            }, fh, indent=4)

args = parse_args()
check_cutoff_type(args)
reads = read_length_file(args.length)

pairs, hits = read_infile(args, reads)
write_file(args, pairs)
write_id_file(args, hits)

