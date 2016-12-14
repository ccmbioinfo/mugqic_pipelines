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

def vardict_varcall(input, output, sample_name, bed_file):

    return Job(
        [input],
        [output],
        [
            ['vardict', 'module_vardict'],
	    ['vardict', 'module_R']
        ],
        command="""\
VarDict -G {reference_fasta} {vardict_options} -th 12 -N {sample_name} \\
  -b {input} -c 1 -S 2 -E 3 -g 4 {bed_file} | teststrandbias.R | var2vcf_valid.pl -N {sample_name} \\
 -E {vardict_options} {output}""".format(
	reference_fasta=config.param('vardict','genome_fasta',type='filepath'),
	vardict_options=config.param("vardict","vardict_other_options"),
	sample_name=sample_name,
	bed_file=bed_file,
        input=" \\\n " + input if input else "",
        output=" \\\n  > " + output if output else ""
        )
    )


# Due to memory issues with running the VarDict on the entire genome, the following method runs VarDict on individual chromosomes
def vardict_varcall_indivChr(input_bam, output, sample_name):

    return Job(
        [input_bam],
        [output],
        [
            ['vardict', 'module_vardict'],
	    ['vardict', 'module_R']
        ],
        command="""\
myfile=$(echo "{input_bam}" | sed 's/.bam//') && \\
for chr in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y; do \\
    VarDict -G {reference_fasta} {vardict_options} -th 12 -N {sample_name} \\
  -b $myfile.REF_$chr.bam -c 1 -S 2 -E 3 -g 4 /hpf/largeprojects/ccmbio/samkh/bioinf_rep/bedfile_test/chr$chr.bed | teststrandbias.R | var2vcf_valid.pl -N {sample_name} \\
 -E {vardict_options} $myfile.vardict.REF_$chr.vcf; done""".format(
	reference_fasta=config.param('vardict','genome_fasta',type='filepath'),
	vardict_options=config.param("vardict","vardict_other_options"),
	sample_name=sample_name,
        input_bam=" \\\n " + input_bam if input_bam else "",
        output=" \\\n  > " + output if output else ""
        )
    )


#-b $myfile.REF_$chr.bam -z -R chr1:125000-58270348 | teststrandbias.R | var2vcf_valid.pl -N {sample_name} \\
