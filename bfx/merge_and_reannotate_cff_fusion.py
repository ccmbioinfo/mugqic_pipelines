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



def merge_and_reannotate_cff_fusion(input_cff, out_dir, annotation_file=None, reference_file=None, ini_section='merge_and_reannotate_cff_fusion'):
	other_options = config.param(ini_section, 'other_options', required=False)
	merged_cff = os.path.join(out_dir, "merged.cff")
	return Job(
		input_cff,
		[os.path.join(out_dir,  merged_cff+".reann")],
		command="""\
cat {cff_files} > {out_dir}/merged.cff &&
/hpf/largeprojects/ccmbio/jiangyue/DIPG_analysis_by_samples/Scripts/reann_cff_fusion.py {merged_cff} {annotation_file} {reference_file} > {merged_cff}.reann && 
/hpf/largeprojects/ccmbio/jiangyue/DIPG_analysis_by_samples/Analysis/generate_common_fusion_stats.py {merged_cff}.reann > {merged_cff}.reann.cluster
""".format(
		cff_files=" \\\n".join(input_cff),
		out_dir=out_dir,
		merged_cff=merged_cff,
		annotation_file=annotation_file if annotation_file else config.param(ini_section, 'annotation_file', type='filepath'),
		reference_file=reference_file if reference_file else config.param(ini_section, 'reference_file', type='filepath')
		),
		removable_files=[]
	)

