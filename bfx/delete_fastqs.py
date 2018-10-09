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



def delete_fastqs(sample, fusion_result_file, ini_section='delete_fastqs'):
	#defuse_result = os.path.join("fusions", "defuse", sample, "results.filtered.tsv")
	#fusionmap_result = os.path.join("fusions", "fusionmap", sample, "02_RNA.FusionReport.txt")
	#ericscript_result = os.path.join("fusions", "ericscript", sample, "fusion.results.filtered.tsv")
	#integrate_result = os.path.join("fusions", "integrate", sample, "breakpoints.tsv")

	out_file=os.path.join("delete_fastqs", "done")

	return Job(
		fusion_result_file,
		[out_file],
		[],
		command="""\
rm -rf {fastq_folder} && rm -f {tophat2_bam} && touch {out_file}""".format(
		fastq_folder=os.path.join("fusions", "gunzip_fastq", sample),
		eric_out=os.path.join("fusions", "ericscript", sample, "out"),
		tophat2_bam=os.path.join("fusions", "tophat2", sample, "*.ba?"),
		out_file=out_file
		),
		removable_files=[]
	)

