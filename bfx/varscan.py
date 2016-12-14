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

def mpileupcns(input, output, sampleNamesFile, other_options=None):

    return Job(
        [input, sampleNamesFile],
        [output],
        [
            ['varscan', 'module_java'],
            ['varscan', 'module_varscan'],
        ],
        command="""\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar $VARSCAN_JAR mpileup2cns {other_options}{input} \\
  --output-vcf 1 \\
  --vcf-sample-list {sampleNames}{output}""".format(
        tmp_dir=config.param('varscan', 'tmp_dir'),
        java_other_options=config.param('varscan', 'java_other_options'),
        ram=config.param('varscan', 'ram'),
        other_options=other_options,
        input=" \\\n " + input if input else "",
        sampleNames=sampleNamesFile,
        output=" \\\n  > " + output if output else ""
        )
    )


def mpileupcns_jacek(input, output, other_options=None):

    return Job(
        [input],
        [output],
        [
            ['varscan', 'module_java'],
            ['varscan', 'module_varscan'],
        ],
        command="""\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar $VARSCAN_JAR mpileup2cns {input}{other_options}{output}""".format(
        tmp_dir=config.param('varscan', 'tmp_dir'),
        java_other_options=config.param('varscan', 'java_other_options'),
        ram=config.param('varscan', 'ram'),
        input="\\\n  " + input if input else "",
        other_options=" \\\n  " + other_options,
        output=" \\\n  > " + output if output else ""
        )
    )


def mpileupcns_jacek_indivChr(input_bam, output, other_options=None):

    return Job(
        [input_bam],
        [output],
        [
            ['varscan', 'module_java'],
            ['varscan', 'module_varscan'],
        ],
        command="""\
myfile=$(echo "{input_bam}" | sed 's/.bam//') && \\
for chr in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y MT; do \\
java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} -jar $VARSCAN_JAR mpileup2cns $myfile.REF_$chr.bam{other_options} $myfile.REF_$chr_varscan.vcf;
done""".format(
        tmp_dir=config.param('varscan', 'tmp_dir'),
        java_other_options=config.param('varscan', 'java_other_options'),
        ram=config.param('varscan', 'ram'),
        input_bam="\\\n  " + input_bam if input_bam else "",
        other_options=" \\\n  " + other_options,
        output=" \\\n  > " + output if output else ""
        )
    )
