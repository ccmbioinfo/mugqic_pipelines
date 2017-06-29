"""
Extract tophits from alignment output (BWA,BLAT,DIAMOND)
"""

import json
from argparse import ArgumentParser

def parse_args():
    arg_parser = ArgumentParser()
    arg_parser.add_argument('--length', help='seqeunce length file')
    arg_parser.add_argument('--infile', help='input file')
    arg_parser.add_argument('--outfile', help='file to write output to')
#    arg_parser.add_argument('--json', help='json to store id for future processing')
    arg_parser.add_argument('--id-file', help='txt to store id for future processing')
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
    """
    Read id-length file and store extracted information in dictionary
    """
    reads = {}

    if('nr_all' in length_file):
        with open(length_file, 'r+') as fh:
            for line in fh:
                line = line.strip()
                # Skip header line or empty line
                if(line.startswith('@') or not line): continue

                data = line.split('\t')

                _id = data[0]

                if('nr_all' in length_file):
                    length = int(data[2])
                else:
                    length = int(data[1])

                reads[_id] = length

    return reads

def read_infile(args, reads):
    """
    Read input file and extract tophits
    """
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

            if(dbID == line1):
                k += 1
            else:
                k = 0

            if(k == 0):
                line1 = dbID
                max_score = float(line[11])

                hitgiID = line[1]

                pident = float(line[2])
                pident = '{:.2f}'.format(pident)

                poverlap = float(line[3])
                poverlap = '{:.2f}'.format(poverlap)

                evalue = line[10]

                if(int(args.cutoff_type) == 0):
                    pairs[dbID] = '\t'.join(str(x) for x in (dbID, hitgiID, pident, poverlap, evalue, max_score))
                    hits[dbID] = hitgiID
                else:
                    if(dbID not in reads):
                        reads[dbID] = 100
                    if(args.diamond):
                        poverlap = 100 * (3 * int(line[3])) / reads[dbID]
                    else:
                        poverlap = 100 * int(line[3]) / reads[dbID]

                    poverlap = '{:.2f}'.format(float(poverlap))

                    if(reads[dbID] >= int(args.cutoff0)):
                        if(max_score >= int(args.cutoff3)):
                            if(dbID in pairs):
                                tmp = pairs[dbID].split('\t')
                                if(tmp[5] < max_score):
                                    pairs[dbID] = '\t'.join(str(x) for x in (dbID, hitgiID, pident, poverlap, evalue, max_score))
                            else:
                                pairs[dbID] = '\t'.join(str(x) for x in (dbID, hitgiID, pident, poverlap, evalue, max_score))
                            hits[dbID] = hitgiID
                        else:
                            max_score += 1
                    else:
                        if((pident >= int(args.cutoff1)) and
                            (poverlap >= int(args.cutoff2))):
                            if(dbID in pairs):
                                tmp = pairs[dbID].split('\t')
                                if(tmp[5] < max_score):
                                    pairs[dbID] = '\t'.join(str(x) for x in (dbID, hitgiID, pident, poverlap, evalue, max_score))
                            else:
                                pairs[dbID] = '\t'.join(str(x) for x in (dbID, hitgiID, pident, poverlap, evalue, max_score))
                            hits[dbID] = hitgiID
                        else:
                            max_score += 1
    return pairs, hits

def write_file(args, pairs):
    """
    Write txt file with all the information
    """
    with open(args.outfile, 'w+') as fh:
        for key in pairs:
            fh.write(pairs[key] + '\n')

def write_id_file(args, hits):
    """
    Write ID and alignment hit to file
    """
    with open(args.id_file, 'w+') as fh:
        fh.write("@id\thit\n")
        for _id in hits:
            fh.write("{_id}\t{hit}\n".format(_id=_id,
                                            hit=hits[_id]))


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

