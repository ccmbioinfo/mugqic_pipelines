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

def fastqc(fastq_files, output_directory, ini_section='fastqc'):
    out_files = [] # out_files for the job parameter
    
    for filename in fastq_files:
        # Assumes the file ends in fastq.gz
        out_name_prefix = os.path.join(output_directory, os.path.basename('.'.join(filename.split('fastq')[0].split('.')[:-1])))
        out_files.extend([out_name_prefix + "_fastqc.html", out_name_prefix + "_fastqc.zip"])

    return Job(
        fastq_files,
        out_files,
        [
            ['fastqc', 'module_fastqc']
        ],
        command="""\
fastqc -t {nb_threads} {in_file_list} -o {out_dir} --nogroup""".format(
            nb_threads=config.param(ini_section, 'threads'),
            in_file_list=' '.join(fastq_files),
            out_dir=output_directory
            )
        )
