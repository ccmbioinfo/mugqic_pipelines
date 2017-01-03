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



def cff_convert(sample, fusion_result_file, sample_type, disease_name, tool, out_dir, ini_section='cff_convertion'):
	#defuse_result = os.path.join("fusions", "defuse", sample, "results.filtered.tsv")
	#fusionmap_result = os.path.join("fusions", "fusionmap", sample, "02_RNA.FusionReport.txt")
	#ericscript_result = os.path.join("fusions", "ericscript", sample, "fusion.results.filtered.tsv")
	#integrate_result = os.path.join("fusions", "integrate", sample, "breakpoints.tsv")


	return Job(
		[fusion_result_file],
		[os.path.join(out_dir, sample+"."+tool+".cff")],
		[],
		command="""\
/hpf/largeprojects/ccmbio/jiangyue/DIPG_analysis_by_samples/Scripts/convert_fusion_results_to_cff.py {sample} {sample_type} {disease_name} {tool} {fusion_result_file} {out_dir}
""".format(
		sample=sample,
		sample_type=sample_type,
		disease_name=disease_name,
		fusion_result_file=fusion_result_file,
		tool=tool,
		out_dir=out_dir
		),
		removable_files=[]
	)

