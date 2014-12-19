#!/usr/bin/env python

# Python Standard Modules
import csv
import logging
import os
import re

# MUGQIC Modules
from sample import *

log = logging.getLogger(__name__)

class Contrast:

    def __init__(self, name):
        self._name = name
        self._controls = []
        self._treatments = []

    @property
    def name(self):
        return self._name

    @property
    def controls(self):
        return self._controls

    @property
    def treatments(self):
        return self._treatments


def parse_new_design_file(design_file, samples):

    log.info("Parse design file " + design_file + " ...")
    design_csv = csv.DictReader(open(design_file, 'rb'), delimiter='\t')

    # Skip first column which is Sample
    contrasts = [Contrast(name) for name in design_csv.fieldnames[1:]]

    for line in design_csv:

        sample_name = line['Sample']
        matching_samples = [sample for sample in samples if sample.name == sample_name]
        if matching_samples:
            # There should be only one matching sample
            sample = matching_samples[0]
        else:
            raise Exception("Error: sample " + sample_name + " in design file " + design_file + " not found in pipeline samples!")

        # Skip first column which is Sample
        for contrast in contrasts:
            sample_contrast_type = line[contrast.name]
            # Empty types are ignored
            if sample_contrast_type:
                if sample_contrast_type == "control":
                    contrast.controls.append(sample)
                elif sample_contrast_type == "treatment":
                    contrast.treatments.append(sample)
                else:
                    raise Exception("Error: invalid value for sample " + sample_name + " and contrast " + contrast.name + " in design file " + design_file + " (should be 'control', 'treatment' or '')!")

    for contrast in contrasts:
        log.info("Contrast " + contrast.name + " (controls: " + str(len(contrast.controls)) + ", treatments: " + str(len(contrast.treatments)) + ") created")
    log.info(str(len(contrasts)) + " contrast" + ("s" if len(contrasts) > 1 else "") + " parsed\n")

    return contrasts

def parse_design_file(design_file, samples):

    design_csv = csv.DictReader(open(design_file, 'rb'), delimiter='\t')

    # Skip first column which is Sample
    contrasts = [Contrast(name) for name in design_csv.fieldnames[1:]]

    for line in design_csv:

        sample_name = line['Sample']
        matching_samples = [sample for sample in samples if sample.name == sample_name]
        if matching_samples:
            # There should be only one matching sample
            sample = matching_samples[0]
        else:
            raise Exception("Error: sample " + sample_name + " in design file " + design_file + " not found in pipeline samples!")

        for contrast in contrasts:
            sample_contrast_type = line[contrast.name]
            # Empty or '0' types are ignored
            if not sample_contrast_type or sample_contrast_type == "0":
                pass
            elif sample_contrast_type == "1":
                contrast.controls.append(sample)
            elif sample_contrast_type == "2":
                contrast.treatments.append(sample)
            else:
                raise Exception("Error: invalid value for sample " + sample_name + " and contrast " + contrast.name + " in design file " + design_file + " (should be '1' for control, '2' for treatment, '0' or '' to be ignored)!")

    for contrast in contrasts:
        log.info("Contrast " + contrast.name + " (controls: " + str(len(contrast.controls)) + ", treatments: " + str(len(contrast.treatments)) + ") created")
    log.info(str(len(contrasts)) + " contrast" + ("s" if len(contrasts) > 1 else "") + " parsed\n")

    return contrasts
