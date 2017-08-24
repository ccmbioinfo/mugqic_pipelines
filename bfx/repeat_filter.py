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



def repeat_filter(cff_file, out_dir, seq_len=None, ref_file=None, config_file=None, ini_section='repeat_filter'):
	other_options = config.param(ini_section, 'other_options', required=False)
	seq_len=seq_len if seq_len else config.param(ini_section, 'seq_len', type='int')
	result_file = os.path.join(".".join([cff_file, "bwafilter", str(seq_len)]))
	return Job(
		[cff_file, ref_file, config_file if config_file else None],
		[result_file],
		[["bwa", "module_bwa"],['cff_convertion', 'module_fusiontools']],
		command="""\
fusion_gene_seq_to_fasta.py {cff_file} {seq_len} > {cff_file}.fa && 
bwa mem {bwa_idx} {cff_file}.fa > {cff_file}.fa.sam &&
filter_fusion_on_bwa_aln.py {cff_file} {cff_file}.fa.sam > {cff_file}.bwafilter.{seq_len} && 
rm {cff_file}.fa && rm {cff_file}.fa.sam""".format(
		other_options=" \\\n  " + other_options if other_options else "",
		cff_file=cff_file,
		seq_len=seq_len,
		bwa_idx=ref_file if ref_file else config.param(ini_section, 'genome_bwa_index', type='filepath'),
		),
		removable_files=[os.path.join(out_dir, cff_file + ".fa"), os.path.join(out_dir, cff_file + ".fa.sam")]
		

	)

