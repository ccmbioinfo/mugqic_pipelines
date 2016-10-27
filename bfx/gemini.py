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

#!/usr/bin/env python

# Python Standard Modules

# MUGQIC Modules
from core.config import *
from core.job import *


def gemini_annotations(variants, gemini_output, tmp_dir, annotations=''):
    """
    Store variants in a gemini database

    :param variants: vcf file to load
    :param gemini_output: db file to load to
    :param tmp_dir: dir to put temporary database chunks
    :param annotations: which annotations were use, either 'VEP', or 'snpEff'
    :return: Job
    """
    # Temporarily remove PYTHONPATH so that we import platform from gemini
    return Job(
        [variants],
        [gemini_output],
        [
            ['gemini_annotations', 'module_gemini'],
            ['gemini_annotations', 'module_htslib'],
            ['gemini_annotations', 'module_tabix']
        ],
        command="""\
PYTHONPATH='' gemini load -v {variants} \\
{options} {annotations} \\
--tempdir {temp} \\
{output}""".format(
        options=config.param('gemini_annotations', 'options'),
        annotations='-t ' + annotations if annotations else '',
        variants=variants,
        output=gemini_output,
        temp=tmp_dir
        )
    )
