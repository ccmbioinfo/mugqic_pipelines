#!/usr/bin/env python

# Python Standard Modules
import csv
import logging
import os
import re

# MUGQIC Modules
from sample import *

log = logging.getLogger(__name__)

class Readset(object):

    def __init__(self, name):
        if re.search("^\w[\w.-]*$", name):
            self._name = name
        else:
            raise Exception("Error: readset name \"" + name +
                "\" is invalid (should match [a-zA-Z0-9_][a-zA-Z0-9_.-]*)!")

    @property
    def name(self):
        return self._name

    @property
    def sample(self):
        return self._sample


class IlluminaReadset(Readset):

    def __init__(self, name, run_type):
        super(IlluminaReadset, self).__init__(name)

        if run_type in ("PAIRED_END", "SINGLE_END"):
            self._run_type = run_type
        else:
            raise Exception("Error: readset run_type \"" + run_type +
                "\" is invalid (should be \"PAIRED_END\" or \"SINGLE_END\")!")

        self.fastq1 = None
        self.fastq2 = None

    @property
    def run_type(self):
        return self._run_type

    @property
    def bam(self):
        if not hasattr(self, "_bam"):
            return None
        else:
            return self._bam

    @property
    def library(self):
        return self._library

    @property
    def run(self):
        return self._run

    @property
    def lane(self):
        return self._lane

    @property
    def adaptor1(self):
        return self._adaptor1

    @property
    def adaptor2(self):
        return self._adaptor2

    @property
    def quality_offset(self):
        return self._quality_offset

    @property
    def beds(self):
        return self._beds

def parse_illumina_readset_file(illumina_readset_file):
    readsets = []
    samples = []

    log.info("Parse Illumina readset file " + illumina_readset_file + " ...")
    readset_csv = csv.DictReader(open(illumina_readset_file, 'rb'), delimiter='\t')
    for line in readset_csv:
        sample_name = line['Sample']
        sample_names = [sample.name for sample in samples]
        if sample_name in sample_names:
            # Sample already exists
            sample = samples[sample_names.index(sample_name)]
        else:
            # Create new sample
            sample = Sample(sample_name)
            samples.append(sample)

        # Create readset and add it to sample
        readset = IlluminaReadset(line['Readset'], line['RunType'])

        # Readset file paths are either absolute or relative to the readset file
        # Convert them to absolute paths
        for format in ("BAM", "FASTQ1", "FASTQ2"):
            if line.get(format, None):
                line[format] = os.path.expandvars(line[format])
                if not os.path.isabs(line[format]):
                    line[format] = os.path.dirname(os.path.abspath(os.path.expandvars(illumina_readset_file))) + os.sep + line[format]
                line[format] = os.path.normpath(line[format])

        readset._bam = line.get('BAM', None)
        readset.fastq1 = line.get('FASTQ1', None)
        readset.fastq2 = line.get('FASTQ2', None)
        readset._library = line.get('Library', None)
        readset._run = line.get('Run', None)
        readset._lane = line.get('Lane', None)
        readset._adaptor1 = line.get('Adaptor1', None)
        readset._adaptor2 = line.get('Adaptor2', None)
        readset._quality_offset = int(line['QualityOffset']) if line.get('QualityOffset', None) else None
        readset._beds = line['BED'].split(";") if line.get('BED', None) else []

        readsets.append(readset)
        sample.add_readset(readset)

    log.info(str(len(readsets)) + " readset" + ("s" if len(readsets) > 1 else "") + " parsed")
    log.info(str(len(samples)) + " sample" + ("s" if len(samples) > 1 else "") + " parsed\n")
    return readsets

class IlluminaRawReadset(IlluminaReadset):

    def __init__(self, name, run_type):
        super(IlluminaRawReadset, self).__init__(name, run_type)

    @property
    def index(self):
        return self._index

    @property
    def reference_species(self):
        return self._reference_species

    @property
    def reference_assembly(self):
        return self._reference_assembly

    @property
    def aligner(self):
        return self._aligner

    @property
    def aligner_reference_file(self):
        return self._aligner_reference_file

    @property
    def reference_file(self):
        return self._reference_file

    @property
    def species(self):
        return self._species

    @property
    def project(self):
        return self._project

    @property
    def library_source(self):
        return self._library_source

    @property
    def operator(self):
        return self._operator

    @property
    def recipe(self):
        return self._recipe

    @property
    def control(self):
        return self._control

    @property
    def description(self):
        return self._description

    @property
    def flow_cell(self):
        return self._flow_cell


def parse_illumina_raw_readset_files(output_dir, run_type, nanuq_readset_file, casava_sheet_file, lane, default_species_genome, genome_root):
    readsets = []
    samples = []

    # Parsing Nanuq readset sheet
    log.info("Parse Nanuq Illumina readset file " + nanuq_readset_file + " ...")
    readset_csv = csv.DictReader(open(nanuq_readset_file, 'rb'), delimiter=',', quotechar='"')
    for line in readset_csv:
        current_lane = line['Region']

        if (int(current_lane) != lane):
            continue

        sample_name = line['Name']

        # Always create a new sample
        sample = Sample(sample_name)
        samples.append(sample)

        # Create readset and add it to sample
        readset = IlluminaRawReadset(line['ProcessingSheetId'], run_type)
        readset._quality_offset = 33
        readset._library = line['Library Barcode']
        readset._library_source = line['Library Source']

        # TODO change aligner to support RNA-seq data
        #if re.search("RNA|cDNA", readset.library_source):
        #    readset._aligner = "star"
        #else:
        #    readset._aligner = "bwa"
        readset._aligner = "bwa"

        readset._run = line['Run']
        readset._lane = current_lane
        readset._adaptor1 = line['Adaptor Read 1 (NOTE: Usage is bound by Illumina Disclaimer found on Nanuq Project Page)']
        readset._adaptor2 = line['Adaptor Read 2 (NOTE: Usage is bound by Illumina Disclaimer found on Nanuq Project Page)']
        if line['BED Files']:
            readset._beds = line['BED Files'].split(";")
        else:
            readset._beds = []

        readsets.append(readset)
        sample.add_readset(readset)


    # Parsing Casava sheet
    log.info("Parsing Casava sample sheet " + casava_sheet_file + " ...")
    casava_csv = csv.DictReader(open(casava_sheet_file, 'rb'), delimiter=',')
    for line in casava_csv:
        if (int(line['Lane']) != lane):
            continue
        processingSheetId = line['SampleID']
        readset = [x for x in readsets if x.name == processingSheetId][0]
        readset._flow_cell = line['FCID']
        readset._species = line['SampleRef']
        readset._index = line['Index']
        readset._description = line['Description']
        readset._control = line['Control']
        readset._recipe = line['Recipe']
        readset._operator = line['Operator']
        readset._project = line['SampleProject']

        fastq_file_pattern = os.path.join(output_dir,
                                          "Unaligned." + readset.lane,
                                          'Project_' + readset.project,
                                          'Sample_' + readset.name,
                                          readset.name + '_' + readset.index + '_L00' + readset.lane + '_R{read_number}_001.fastq.gz')
        readset.fastq1 = fastq_file_pattern.format(read_number = 1)
        readset.fastq2 = fastq_file_pattern.format(read_number = 2) if readset.run_type == "PAIRED_END" else None;


    # Searching for a matching reference for the specified species
    for readset in readsets:
        # Find if any reference_assembly or reference_species for the specied species
        for genome in default_species_genome.split('~'):
            values = genome.split(':')
            if (re.match(values[0], readset.species, re.IGNORECASE)) :
                aligner_reference_file = ""
                if (readset.aligner == "bwa"):
                    aligner_reference_file = os.path.join(genome_root, values[1] + "." + values[2],
                                                          "genome",
                                                          "bwa_index",
                                                          values[1] + "." + values[2] + ".fa")
                reference_file = os.path.join(genome_root, values[1] + "." + values[2],
                                              "genome",
                                              values[1] + "." + values[2] + ".fa")
                if reference_file and os.path.isfile(reference_file):
                    if aligner_reference_file and os.path.isfile(aligner_reference_file):
                        readset._aligner_reference_file = aligner_reference_file
                        readset._reference_file = reference_file
                        readset._reference_species = values[1]
                        readset._reference_assembly = values[2]
                        readset._bam = os.path.join(output_dir,
                                                    "Aligned." + readset.lane,
                                                    'alignment',
                                                    readset.sample.name,
                                                    'run' + readset.run + "_" + readset.lane,
                                                    readset.sample.name + readset.library)
                    else:
                        log.warning("Unable to access the aligner reference file: '" + aligner_reference_file + "'")
                else:
                    log.warning("Unable to access the reference file: '" + reference_file + "'")

        if (readset.bam is None):
            log.info("Skipping alignment for the species: " + readset.species)


    log.info(str(len(readsets)) + " readset" + ("s" if len(readsets) > 1 else "") + " parsed")
    log.info(str(len(samples)) + " sample" + ("s" if len(samples) > 1 else "") + " parsed\n")
    return readsets


class PacBioReadset(Readset):

    @property
    def run(self):
        return self._run

    @property
    def smartcell(self):
        return self._smartcell

    @property
    def protocol(self):
        return self._protocol

    @property
    def nb_base_pairs(self):
        return self._nb_base_pairs

    @property
    def estimated_genome_size(self):
        if self._estimated_genome_size:
            return self._estimated_genome_size
        else:
            raise Exception("Error: readset \"" + self.name + "\" estimated_genome_size is not defined!")

    @property
    def bas_files(self):
        return self._bas_files

    @property
    def bax_files(self):
        return self._bax_files

def parse_pacbio_readset_file(pacbio_readset_file):
    readsets = []
    samples = []

    log.info("Parse PacBio readset file " + pacbio_readset_file + " ...")
    readset_csv = csv.DictReader(open(pacbio_readset_file, 'rb'), delimiter='\t')
    for line in readset_csv:
        sample_name = line['Sample']
        sample_names = [sample.name for sample in samples]
        if sample_name in sample_names:
            # Sample already exists
            sample = samples[sample_names.index(sample_name)]
        else:
            # Create new sample
            sample = Sample(sample_name)
            samples.append(sample)

        # Create readset and add it to sample
        readset = PacBioReadset(line['Readset'])

        # Readset file paths are either absolute or relative to the readset file
        # Convert them to absolute paths
        for format in ("BAS", "BAX"):
            if line.get(format, None):
                abs_files = []
                for file in line[format].split(","):
                    file = os.path.expandvars(file)
                    if not os.path.isabs(file):
                        file = os.path.dirname(os.path.abspath(os.path.expandvars(pacbio_readset_file))) + os.sep + file
                    abs_files.append(os.path.normpath(file))
                line[format] = ",".join(abs_files)

        readset._run = line.get('Run', None)
        readset._smartcell = line.get('Smartcell', None)
        readset._protocol = line.get('Protocol', None)
        readset._nb_base_pairs = int(line['NbBasePairs']) if line.get('NbBasePairs', None) else None
        readset._estimated_genome_size = int(line['EstimatedGenomeSize']) if line.get('EstimatedGenomeSize', None) else None
        readset._bas_files = line['BAS'].split(",") if line.get('BAS', None) else []
        readset._bax_files = line['BAX'].split(",") if line.get('BAX', None) else []

        readsets.append(readset)
        sample.add_readset(readset)

    log.info(str(len(readsets)) + " readset" + ("s" if len(readsets) > 1 else "") + " parsed")
    log.info(str(len(samples)) + " sample" + ("s" if len(samples) > 1 else "") + " parsed\n")
    return readsets
