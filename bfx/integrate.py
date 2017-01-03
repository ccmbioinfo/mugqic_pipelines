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



def integrate(accepted_bam, unmapped_bam, out_dir, ini_section='integrate'):
	other_options = config.param(ini_section, 'other_options', required=False)
	result_file = os.path.join(out_dir, "breakpoints.tsv")
	return Job(
		[accepted_bam, unmapped_bam],
		[result_file],
		[["bwa", "module_bwa"]],
		command="""\
/hpf/largeprojects/ccmbio/jiangyue/tools/INTEGRATE_0_2_0/Integrate/INTEGRATE-build/bin/Integrate fusion {other_options} /hpf/largeprojects/ccmbio/jiangyue/hg19_decoy/human_g1k_v37_decoy.fasta /hpf/largeprojects/ccmbio/jiangyue/tools/INTEGRATE_0_2_0/Integrate/annot.ucsc.txt /hpf/largeprojects/ccmbio/jiangyue/hg19_decoy/integrate_index {accepted_bam} {unmapped_bam}""".format(
		other_options=" \\\n  " + other_options if other_options else "",
		accepted_bam=accepted_bam,
		unmapped_bam=unmapped_bam
		),
		removable_files=[]
	)

