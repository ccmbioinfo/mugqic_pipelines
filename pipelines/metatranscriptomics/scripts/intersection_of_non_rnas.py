#!/usr/bin/env python

from Bio import SeqIO


def pre_slash(id):
    return id.split('/')[0]


cow = dict()
for i in (1, 2):
    original_fastq = 'remove_rrna/cow{}_qual_unique_n_rRNA.fastq'.format(i)
    print('Getting reads from ' + original_fastq)
    cow[i] = set(SeqIO.parse(original_fastq, 'fastq'))

print('Taking intersection of ids')
ids_in_both = {pre_slash(e.id) for e in cow[1]} & {pre_slash(e.id) for e in cow[2]}

cow_intersection = dict()
for i in (1, 2):
    print('Filtering cow' + str(i))
    cow_intersection[i] = sorted([e for e in cow[i] if pre_slash(e.id) in ids_in_both], key = lambda e: pre_slash(e.id))

    out_file = 'remove_host_reads/cow{}_matching_ids.fastq'.format(i)
    print('Writing sequences to ' + out_file)
    with open(out_file, 'w+') as f:
        SeqIO.write(cow_intersection[i], f, 'fastq')
