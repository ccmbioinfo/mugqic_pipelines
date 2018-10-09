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

def cff_convert(sample, fusion_result_file, sample_info_file, tool, out_dir, ini_section='cff_conversion'):

    return Job(
        [fusion_result_file],
        [os.path.join(out_dir, sample+"."+tool+".cff")],
        [['cff_conversion', 'module_fusiontools']],
        command="""\
convert_fusion_results_to_cff.py \\
  {sample} \\
  {sample_info_file} \\
  {tool} \\
  {fusion_result_file} \\
  {out_dir}""".format(
        sample=sample,
        sample_info_file=sample_info_file,
        fusion_result_file=fusion_result_file,
        tool=tool,
        out_dir=out_dir
        ),
        removable_files=[]
    )

