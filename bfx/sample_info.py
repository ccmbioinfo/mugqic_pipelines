#!/usr/bin/env python

################################################################################
# Copyright (C) 2014, 2015 GenAP, McGill University and Genome Quebec Innovation Centre
#
# This file is part of MUGQIC Pipelines.
#
# MUGQIC Pipelines is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# MUGQIC Pipelines is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with MUGQIC Pipelines.  If not, see <http://www.gnu.org/licenses/>.
################################################################################

# Python Standard Modules
import logging
import os

log = logging.getLogger(__name__)


# TODO: Clean this up...
def parse_sample_info(info_file):

    log.info("Parse sample info file " + info_file + " ...")
    sample_info = {}

    with open(sample_info_file, "r") as f:
        for line in f.readlines():
            line = line.strip()
            # Ignore comments
            if not line.startswith("#"):
                split_sample = line.split("\t")
                # Make sure the format is correct
                if len(split_sample) > 1:
                    param = split_sample[1].split("=")
                    if split_sample[0] in sample_info.keys():
                        sample_info[split_sample[0]][param[0]] = param[1] if len(param) > 1 else 1 # If not equal to anything, set to true b/c called
                    else:
                        sample_info[split_sample[0]] = {param[0]: param[1]} if len(param) > 1 else {param[0]: 1}

    return sample_info
