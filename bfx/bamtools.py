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

# MUGQIC Modules
from core.config import *
from core.job import *

def split(input_bam, output, other_options=None):

    return Job(
        [input_bam],
        [output],
        [
            ['bamtools', 'module_bamtools']
        ],
        command="""\
bamtools split -in {input_bam}  \\
  -reference && touch {output}""".format(
        tmp_dir=config.param('bamtools', 'tmp_dir'),
        ram=config.param('bamtools', 'ram'),
        other_options=other_options,
        input_bam=" \\\n " + input_bam if input_bam else "",
        output=" \\\n  " + output if output else ""
        )
    )


def index(input_bam, bamtools_split_done, output, other_options=None):

    return Job(
        [input_bam, bamtools_split_done],
        [output],
        [
            ['bamtools', 'module_bamtools']
        ],
        command="""\
for chr in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y MT; do
    myfile=$(echo "{input_bam}" | sed 's/.bam//') && \\
    bamtools index -in $myfile.REF_$chr.bam;  \\
done && touch {output}""".format(
        tmp_dir=config.param('bamtools', 'tmp_dir'),
        ram=config.param('bamtools', 'ram'),
        other_options=other_options,
        input_bam=" \\\n " + input_bam if input_bam else "",
        output=" \\\n  " + output if output else ""
        )
    )


