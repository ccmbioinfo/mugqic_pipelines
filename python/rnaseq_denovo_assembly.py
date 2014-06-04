#!/usr/bin/env python

# Python Standard Modules
import argparse
import collections

# MUGQIC Modules
from pipeline import *
from readset import *
from trimmomatic import *


class RnaSeqDeNovoAssembly(Pipeline):

    @property
    def readsets(self):
        return self._readsets

    @property
    def samples(self):
        return self._samples

    def trim(self, readset):
        return trimmomatic(readset.name, readset.name + ".trim")

    def normalization(self, readset):
        return normalize(readset.name + ".trim", readset.name + ".trim.normalized")

    def align(self):
        return trinity([readset.name + ".trim.normalized" for readset in self.readsets])

    def blast(self):
        return blastx("Trinity.fasta")

    def abundance(self, sample):
        return rsem("Trinity.fasta", sample.name)

    def annotate(self):
        return trinotate("Trinity.fasta")

    def deliverable(self):
        return nozzle(["Trinity.fasta", "Trinity_stats.csv", "trinotate.tsv"] + ["rsem_" + sample.name + ".fpkm" for sample in self.samples], "report")

    @property
    def step_dict_map(self):
        return [
            {"name": self.trim, "loop": self.readsets},
            {"name": self.normalization, "loop": self.readsets},
            {"name": self.align},
            {"name": self.blast},
            {"name": self.abundance, "loop": self.samples},
            {"name": self.annotate},
            {"name": self.deliverable}
        ]

    def __init__(self):
        # Initialize attributes to avoid AttributeError in default_argparser(self.step_dict_map)
        self._readsets = []
        self._samples = []

        argparser = default_argparser(self.step_dict_map)
        # Add pipeline specific arguments
        argparser.add_argument("-r", "--readsets", help="readset file", type=file, required=True)
        argparser.add_argument("-d", "--design", help="design file", type=file)
        args = argparser.parse_args()

        self._readsets = Readset.parse_readset_file(args.readsets.name)
        # Retrieve unique samples from their readsets, removing duplicates
        self._samples = list(collections.OrderedDict.fromkeys([readset.sample for readset in self._readsets]))

        Pipeline.__init__(self, args)
        
RnaSeqDeNovoAssembly().show()