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

def ericscript(in1fastq, in2fastq, out_dir, config_file=None, ini_section='ericscript'):

    other_options = config.param(ini_section, 'other_options', required=False)
    result_file = os.path.join(out_dir, "fusion.results.filtered.tsv")
    ericscript_path = "/hpf/largeprojects/ccmbio/mapostolides/gene_fusion/modules/ericscript-0.5.4/"
    database_path = "/hpf/largeprojects/ccmbio/jiangyue/database/ericscript/ericscript_db_homosapiens_ensembl73"
    #database_path = "/hpf/largeprojects/ccmbio/mapostolides/database/ericscript_db"
    return Job(
        [in1fastq, in2fastq, config_file if config_file else None],
        [result_file],
        [
            ["ericscript", "module_bedtools"],
            ["ericscript", "module_blat"],
            ["ericscript", "module_samtools"],
            ["ericscript", "module_R_3_1_0"],
            ["ericscript", "module_bwa"],
            ["ericscript", "module_seqtk"]

        ],
        command="""\
export PATH=$PATH:{ericscript_path} && \\
ericscript.nodbcheck.pl \\
  {other_options} \\
  -db {database_path} \\
  -name "fusion" \\
  -o {out_dir} \\
  {in1fastq} \\
  {in2fastq} && \\
rm -rf {out_dir}/aln && rm -rf {out_dir}/out""".format(
        ericscript_path=ericscript_path,
        other_options=" \\\n  " + other_options if other_options else "",
        database_path=database_path,
        in1fastq=in1fastq,
        in2fastq=in2fastq,
        out_dir=out_dir
        ),
        removable_files=[]
    )

