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



def merge_and_summary(sum_files, out_dir, ini_section='merge_and_summary'):
	other_options = config.param(ini_section, 'other_options', required=False)
	return Job(
		sum_files,
		[os.path.join(out_dir, "merged.summary.stat")],
		command="""\
cat {sum_files} |grep ^SUM> {out_dir}/merged.summary &&
/hpf/largeprojects/ccmbio/jiangyue/DIPG_analysis_by_samples/Scripts/generate_summary_stats_meregd.py {out_dir}/merged.summary > {out_dir}/merged.summary.stat""".format(
		sum_files=" \\\n".join(sum_files),
		out_dir=out_dir,
		),
		removable_files=[]
	)

