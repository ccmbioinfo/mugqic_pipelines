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
import os
import zipfile
# MUGQIC Modules
from core.config import *
from core.job import *

def gunzip_fastq(infastqgz, out_dir, ini_section='gunzip_fastq'):
    other_options = config.param(ini_section, 'other_options', required=False)
    
    if infastqgz.endswith(".gz"):
        outfastq = os.path.join(out_dir, os.path.basename(os.path.splitext(infastqgz)[0]))
        return Job(
            [infastqgz],
            [outfastq],
            [],
            command="""\
zcat {input} > {output}""".format(
            input=infastqgz,    
            output=outfastq
            )
        )
    else:
        outfastq = os.path.join(out_dir, os.path.basename(infastqgz))
        return Job(
            [infastqgz],
            [outfastq],
            [],
            command="""\
    ln -fs {input}  {output}""".format(
            input=infastqgz,    
            output=outfastq
            )
        )
