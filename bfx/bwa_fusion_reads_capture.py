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
def bwa_fusion_reads_capture(in1fastq, in2fastq, ref, out_bam, read_group=None, ini_section='bwa_fusion_reads_capture'):
	other_options = config.param(ini_section, 'other_options', required=False)

	return Job(
		[in1fastq, in2fastq, ref],
		[out_bam],
		[["bwa_mem", "module_bwa"], ["samtools", "module_samtools"]],
		command="""\
bwa mem {other_options}{read_group} \\
  {ref} \\
  {in1fastq}{in2fastq}|awk 'and($2,0x04)==0 || and($2,0x08)==0'|samtools view -b - > {out_bam} && 
test {pipestatus} -eq 0 &&
samtools sort {out_bam} -o {out_bam} && 
samtools index {out_bam}""".format(
		other_options=" \\\n  " + other_options if other_options else "",
		read_group=" \\\n  -R " + read_group if read_group else "",
		ref=ref,
		in1fastq=in1fastq,
		in2fastq=" \\\n  " + in2fastq,
		out_bam=" \\\n " + out_bam,
		pipestatus="${PIPESTATUS}"
		),
		removable_files=[]
	)
