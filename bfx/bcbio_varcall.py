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

def bcbio_var_recall(input1, input2, input3, input4, input5, output, other_options=None):

    return Job(
        [input1, input2, input3, input4, input5],
        [output],
        [
            ['bcbio_var_recall', 'module_htslib'],
            ['bcbio_var_recall', 'module_bcftools'],
            ['bcbio_var_recall', 'module_tabix']
        ],
        command="""\
{script_path}/scripts/bcbio-variation-recall ensemble --numpass 2 --names GATK_HC,SamtoolsMpileup,Varscan,Vardict,Freebayes \\
  {output} {reference_fasta} {input1}{input2}{input3}{input4}{input5} \\
  """.format(
	script_path=config.param("DEFAULT", "forge_location"),
        reference_fasta=config.param('bcbio_var_recall','genome_fasta',type='filepath'),
        ram=config.param('bcbio_var_recall', 'ram'),
        other_options=other_options,
        input1=" \\\n " + input1 if input1 else "",
        input2=" \\\n " + input2 if input2 else "",
        input3=" \\\n " + input3 if input3 else "",
        input4=" \\\n " + input4 if input4 else "",
        input5=" \\\n " + input5 if input5 else "",
        output=" \\\n  " + output if output else ""
        )
    )


