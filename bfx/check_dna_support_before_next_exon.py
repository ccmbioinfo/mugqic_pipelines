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



def check_dna_support_before_next_exon(input_reann, bam_list_file, tmp_dir, annotation_file=None, ini_section='check_dna_support_before_next_exon'):
	other_options = config.param(ini_section, 'other_options', required=False)
	output_file = input_reann + ".dnasupp"
	return Job(
		[input_reann],
		[output_file],
		[["check_dna_support_before_next_exon", "module_fusiontools"]],
		command="""\
get_fusion_dna_supp_before_end_of_gene.py {reann_file} {bam_list_file} {annotation_file} {tmp_dir} > {output_file}""".format(
		reann_file=input_reann,
		output_file=output_file,
		bam_list_file=bam_list_file,
		tmp_dir=tmp_dir,
		annotation_file=annotation_file if annotation_file else config.param(ini_section, 'annotation_file', type='filepath'),
		),
		removable_files=[os.path.join(tmp_dir, "dicordant.bam"), os.path.join(tmp_dir, "dicordant.bam.bai")]
	)

