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

def tophat2(in1fastq, in2fastq, out_dir, ini_section='tophat2'):

    other_options = config.param(ini_section, 'other_options', required=False)
    accepted_bam = os.path.join(out_dir, "accepted_hits.bam")
    unmapped_bam = os.path.join(out_dir, "unmapped.bam")

    return Job(
        [in1fastq, in2fastq],
        [os.path.join(out_dir, "accepted_hits.bam"), os.path.join(out_dir, "unmapped.bam")],
        [
            ["tophat", "module_tophat"],
            ["tophat", "module_bowtie2"],
            ["tophat", "module_samtools"]
        ],
        command="""\
/hpf/tools/centos6/tophat/2.0.13/tophat2 \\
  {other_options} \\
  -o {out_dir} \\
/hpf/largeprojects/ccmbio/jiangyue/hg19_decoy/bowtie2_index/human_g1k_v37_decoy.fasta \\
  {in1fastq} \\
  {in2fastq} && \\
samtools index {accepted_bam}""".format(
        other_options=" \\\n  " + other_options if other_options else "",
        in1fastq=in1fastq,
        in2fastq=in2fastq,
        out_dir=out_dir,
        accepted_bam=accepted_bam 
        ),
        removable_files=[]
    )

