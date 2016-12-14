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

def freebayes_varcall(input_bam, bamtools_done_file, output=None, other_options=None):

    return Job(
        [input_bam, bamtools_done_file],
        [output],
        [
            ['freebayes', 'module_freebayes']
        ],
        command="""\
myfile=$(echo "{input_bam}" | sed 's/.bam//') && \\
for chr in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y MT; do \\
    freebayes -f {reference_fasta} {freebayes_options} \\
    $myfile.REF_$chr.bam > $myfile.REF_$chr.vcf;  \\
done""".format(
        tmp_dir=config.param('freebayes', 'tmp_dir'),
	reference_fasta=config.param('freebayes','genome_fasta',type='filepath'),
	freebayes_options=config.param("freebayes","freebayes_other_options"),
        ram=config.param('freebayes', 'ram'),
        other_options=other_options,
        input_bam=input_bam,
        bamtools_done_file=bamtools_done_file,
        output=" \\\n  " + output if output else ""
        )
    )



def bcftools_CompressConcat_ChrsVcfs(input_bam, bamtools_done_file, output, other_options=None):

    return Job(
        [input_bam, bamtools_done_file],
        [output],
        [
            ['bcftools', 'module_bcftools'],
	    ['htslib', 'module_htslib'],
	    ['freebayes', 'module_tabix']
        ],
        command="""\
myfile=$(echo "{input_bam}" | sed 's/.bam//') && \\
for chr in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y MT; do \\
    bgzip $myfile.REF_$chr.vcf $myfile.REF_$chr.vcf.gz && \\
    tabix $myfile.REF_$chr.vcf.gz && \\
    vcflist=$vcflist"$myfile.REF_$chr.vcf.gz "; \\
done && \\
bcftools concat $vcflist -o {output} \\
""".format(
        tmp_dir=config.param('bcftools', 'tmp_dir'),
        other_options=other_options,
        input_bam=input_bam,
        bamtools_done_file=bamtools_done_file,
        output=" \\\n  " + output if output else ""
        )
    )


