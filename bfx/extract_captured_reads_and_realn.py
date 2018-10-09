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
def extract_captured_reads_and_realn(captured_bam, realigned_bam, hg_and_trans_junc_ref=None, ini_section='extract_captured_reads_and_realn'):
	other_options = config.param(ini_section, 'other_options', required=False)
	return Job(
		[captured_bam],
		[realigned_bam],
		[["extract_captured_reads_and_realn", "module_fusiontools"], ['bwa', 'module_bwa']],
		command="""\
extract_fusion_reads_DIPG_merged.py {captured_bam} && \\
bwa mem {other_options} \\
  {hg_and_trans_junc_ref} \\
  {captured_bam}.fa > {realigned_bam}""".format(
		other_options=" \\\n  " + other_options if other_options else "",
		captured_bam=captured_bam,
		hg_and_trans_junc_ref=hg_and_trans_junc_ref if hg_and_trans_junc_ref else config.param(ini_section, 'hg_and_trans_junc_ref', type='filepath'),
		realigned_bam=" \\\n " + realigned_bam
		),
		removable_files=[]
	)
