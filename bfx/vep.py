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
from core.config import config
from core.job import Job


def annotate(vcf, out_file):
    config_section = 'gemini_annotations'
    return Job(
        [vcf],
        [out_file],
        [
            [config_section, 'module_perl'],
            [config_section, 'module_tabix']
        ],
        command="""\
perl {vep_location} -i {input} --assembly {assembly} --stats_file {out_file}.summary.html {options} | grep -v -- "- INFO: Disabling" | bgzip -c > {out_file}
""".format(input=vcf,
           vep_location=config.param(config_section, 'vep_location'),
           assembly=config.param(config_section, 'assembly'),
           options=config.param(config_section, 'vep_options'),
           out_file=out_file)
    )


def is_vep_requested():
    return 'vep' in config.param('gemini_annotations', 'annotations').lower()
