#!/hpf/tools/centos6/mugqic-pipelines/source/resource/software/python/Python-2.7.11/bin/python
import re
import os
from Bio import SeqIO


def get_ids(file):
    ids = set()
    with open(file) as f:
        for line in f:
            if line.startswith('S'):
                ids.add(re.match('^(\S+\t){8}(\S+)', line).group(2))
    return ids


for i in (1, 2):
    unique_ids = get_ids('remove_duplicates/cow{}_qual_all_unique.uc'.format(i))
    unique_records = [recrd for recrd in SeqIO.parse('flash/cow{}_qual_all.fastq'.format(i), 'fastq') if recrd.id in unique_ids]

    with open('remove_duplicates/cow{}_qual_all_unique_IDs.txt'.format(i), 'w+') as f:
        for record in unique_records:
           f.write('{id}\t{length}\n'.format(id=id, length=len(record.seq)))

    with open('remove_duplicates/cow{}_qual_all_unique_updated_ids.fasta'.format(i), 'w+') as out:
        SeqIO.write(unique_records, out, 'fasta')
