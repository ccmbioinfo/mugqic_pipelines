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

def build_fusion_genes_ref(cff_file, out_dir, annotation_file=None, reference_file=None, ini_section='build_fusion_and_head_gene_ref'):
	out_file_name=os.path.basename(cff_file) + ".genes.fa"
	return Job(
		[cff_file],
		[os.path.join(out_dir, out_file_name)],
		[['build_fusion_and_head_gene_ref', 'module_fusiontools']],
		command="""\
build_fusion_genes_ref.py {cff_file} {annotation_file} {reference_file} > \\
{out_dir}/{out_file_name} &&
bwa index {out_dir}/{out_file_name}
""".format(
		cff_file=cff_file,
		out_dir=out_dir,
		out_file_name=out_file_name,
		annotation_file=annotation_file if annotation_file else config.param(ini_section, 'annotation_file', type='filepath'),
		reference_file=reference_file if reference_file else config.param(ini_section, 'reference_file', type='filepath')

		)
	)

def build_fusion_and_head_gene_ref(cff_file, out_dir, annotation_file=None, reference_file=None, ini_section='build_fusion_and_head_gene_ref'):
	out_file_name=os.path.basename(cff_file) + ".fa"
	return Job(
		[cff_file],
		[os.path.join(out_dir, out_file_name)],
		[['build_fusion_and_head_gene_ref', 'module_fusiontools'], ['bwa', 'module_bwa']],
		command="""\
build_fusion_and_head_transcript_ref.py {cff_file} {annotation_file} {reference_file} > \\
{out_dir}/{out_file_name} &&
bwa index {out_dir}/{out_file_name}
""".format(
		cff_file=cff_file,
		out_dir=out_dir,
		out_file_name=out_file_name,
		annotation_file=annotation_file if annotation_file else config.param(ini_section, 'annotation_file', type='filepath'),
		reference_file=reference_file if reference_file else config.param(ini_section, 'reference_file', type='filepath')

		)
	)
