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



def defuse(in1fastq, in2fastq, out_dir, config_file=None, ini_section='defuse'):
	other_options = config.param(ini_section, 'other_options', required=False)
	result_file = os.path.join(out_dir, "results.filtered.tsv")
	return Job(
		[in1fastq, in2fastq, config_file if config_file else None],
		[result_file],
		[["defuse", "module_defuse"]],
		command="""\
/hpf/tools/centos6/defuse/0.6.2/scripts/defuse.pl {other_options} -c {config_file}{in1fastq}{in2fastq} -o {out_dir} &&
ls -d {out_dir}/*|grep -v result|xargs rm -rf""".format(
		other_options=" \\\n  " + other_options if other_options else "",
		config_file=" \\\n  " + config_file if config_file else config.param(ini_section, 'defuse_config', required=True),
		in1fastq=" \\\n -1 " + in1fastq,
		in2fastq=" \\\n -2 " + in2fastq,
		out_dir=" \\\n  " + out_dir
		),
		removable_files=[out_dir + "/reads.fqi", out_dir + "/reads.names", out_dir + "/reads.?.fastq"]
		

	)

