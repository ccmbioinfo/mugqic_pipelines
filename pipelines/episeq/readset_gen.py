#!/bin/env python
"""
Info goes here
"""

import csv
import glob
import os
import sys


def parse_manifest(manifest_file):
    # type: (str) -> list(list(str))
    """
    This function will attempt to find the SRA sample number, run number, and the library format.
    This will be sufficient to generate the readset file. The use of DictReader and list comprehension
    allows the function to appear minimal.

    :param manifest_file: The filepath to the NCBI's SRA Run Selector table
    :type manifest_file: str
    :return: Extracted required data from each line
    :rtype: list[list]
    """
    with open(manifest_file) as manifest:
        reader = csv.DictReader(manifest, delimiter='\t', quoting=csv.QUOTE_NONE,
                                skipinitialspace=True, strict=True)
        entries = [[line['SRA_Sample_s'], line['Run_s'],
                    line['LibraryLayout_s'], line['SRA_Study_s']] for line in reader]
    return entries


def generate_readset(entries, readset_file='./episeq.readset', data_root='.'):
    # type: (list, str) -> list
    """
    Writes the readset file, while attempting to find the required fastq/bam files for every run.
    The search algorithm is naive. It uses glob patterns to find any file starting with the run id.
    Bam files are preffered over any number of fastq files. In the event that more than one file are
    eligible, the first file in the list is taken. Unless otherwise specified, the current working
    directory is used as the search location.

    :param data_root: The directory which we use to look for fastq/bam files.
    :type data_root: str
    :param entries: A set of sample information that is parsed from the SRA run table.
    :type entries: list[list]
    :param readset_file: The output readset file to write to.
    :type readset_file: str
    """
    with open(readset_file, 'w') as readsets:
        fieldnames = ['Sample', 'Readset', 'RunType', 'FASTQ1', 'FASTQ2']
        writer = csv.DictWriter(readsets, fieldnames=fieldnames, delimiter='\t', quoting=csv.QUOTE_NONE)
        # writer.writeheader()
        readsets.write('\t'.join(fieldnames) + '\n')

        uniq_study = []

        for entry in entries:
            basename = os.path.join(data_root, entry[1], entry[1])
            found_bam = glob.glob(basename + '*.bam*')
            found_files = glob.glob(basename + '*.fastq*')

            try:
                if found_bam:  # Prefer bam files
                    fastq1 = found_bam[0]
                    fastq2 = ''
                else:  # Handle up to 2 fastq files
                    if len(found_files) == 0:
                        raise IOError('No .fastq* files were found for run: ' + entry[1])
                    elif len(found_files) == 1:
                        fastq1 = found_files[0]
                        fastq2 = ''
                    else:  # Make sure I get the right order... just in case
                        fastq1 = glob.glob(basename + '*_1.fastq*')[0]
                        fastq2 = glob.glob(basename + '*_2.fastq*')[0]
            except IOError:
                print ('%s will not be included in this run because its files are missing.' % entry[1])
                continue  # Don't write anything about missing samples

            writer.writerow({'Sample': entry[0], 'Readset': entry[1], 'RunType': entry[2] + "_END",
                             'FASTQ1': fastq1, 'FASTQ2': fastq2})

            if entry[0] not in [run[0] for run in uniq_study]:
                uniq_study.append([entry[0], entry[3]])

    return uniq_study


def generate_design(study_group, design_file='./episet.design'):
    """
    This function creates/overwrite a design file as an input for the pipeline. It uses information
    from generate_readset() to determine what samples are remaining in the list. Since every
    analysis has a different focus, this function tries to organize samples by study_id and run_number.

    The function only generates one Contrast column and sets each sample to 0. The point is to have
    a specific list of samples that are available for editing.

    :param study_group: A list of tuples. Each tuple must be in the form of (sample_name, Study_id)
    :type study_group: list(list(str))
    :param design_file: The output path and name for the design file. Defaults to CWD.
    :type design_file: str
    """
    with open(design_file, 'w') as design:
        fieldnames = ['Sample', 'Contrast']
        writer = csv.DictWriter(design, fieldnames=fieldnames, restval=0,
                                delimiter='\t', quoting=csv.QUOTE_NONE)
        # writer.writeheader()
        design.write('\t'.join(fieldnames) + '\n')

        study_group = sorted(study_group, key=lambda entry: (entry[1], entry[0]))
        for sample in study_group:
            writer.writerow({'Sample': sample[0], 'Contrast': 0})


if __name__ == '__main__':
    """
    Command line input.
    :arg 1: (Required) A filepath to the SraRunTable.txt or some other similar format.
    :arg 2: (Optional) The output path for the readset file. Default: ./episeq.readset
    :arg 3: (Optional) The output path for the design file. Default : ./episeq.design
    :arg 4: (Optional) The path to the root data directory, where the seq files are stored. Default: '.'
    """
    data = parse_manifest(sys.argv[1])
    study = generate_readset(data, sys.argv[2], sys.argv[4])
    generate_design(study, sys.argv[3])
