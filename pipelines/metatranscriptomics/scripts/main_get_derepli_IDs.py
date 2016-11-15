#!/hpf/tools/centos6/mugqic-pipelines/source/resource/software/python/Python-2.7.11/bin/python
import re
import os
from Bio import SeqIO


def get_ids(file):
    id_to_length = dict()
    with open(file) as f:
        for line in f:
            if line.startswith('S'):
                id_to_length[re.match('^(\S+\t){8}(\S+)', line).group(2)] = None
    return id_to_length


for i in (1, 2):
    id_to_length = get_ids('cow{}_qual_all_unique.uc'.format(i))

    for record in SeqIO.parse('cow{}_qual_all.fastq'.format(i), 'fastq'):
        if record.id in id_to_length:
            id_to_length[record.id] = len(record.seq)

    os.remove('cow{}_qual_all_unique_IDs.txt'.format(i))
    with open('cow{}_qual_all_unique_IDs.txt'.format(i), 'w+') as f:
        for id, length in id_to_length.items():
            f.write('{id}\t{length}\n'.format(id=id, length=length))

    os.remove('cow{}_qual_all_unique.fastq'.format(i))
    with open('cow{}_qual_all_unique.fastq'.format(i), 'w+') as out:
        for record in SeqIO.parse('cow{}_qual_all.fastq'.format(i), 'fastq'):
            if record.id in id_to_length:
                out.write(record.format('fastq'))
