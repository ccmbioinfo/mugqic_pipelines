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

# MUGQIC Modules
from core.config import *
from core.job import *

def annotate_genes(input, output, splicing_threshold, ini_section="annotate_genes"):

    annovar_db_dir = config.param(ini_section, 'annovar_db_dir', type='dirpath')
    other_options = config.param(ini_section, 'other_options', required=False)

    return Job(
        [input],
        [output + ".exonic_variant_function", output + ".variant_function"],
        [
            ['annotate_genes', 'module_annovar'],
            ['annotate_genes', 'module_perl']
        ],
        command="""\
annotate_variation.pl \\
  --buildver {buildver} \\
  --geneanno \\
  --exonsort \\
  --splicing_threshold {splicing_threshold} \\
  --outfile {output}{other_options} \\
  {input} \\
  {annovar_db_dir}""".format(
            buildver=config.param(ini_section, 'assembly'),
            input=input,
            output=output,
            splicing_threshold=splicing_threshold,
            other_options=" \\\n  " + other_options if other_options else "",
            annovar_db_dir=annovar_db_dir
        ))


def combine_annovar_files(vcf_file, exonic_vars, var_fn, var_fn_extended, table_file, output, ini_section="combine_annovar_files"):
    other_options = config.param(ini_section, "other_options", required=False)

    return Job(
        [vcf_file, exonic_vars, var_fn, var_fn_extended, table_file],
        [output],
        [
            ['combine_annovar_files', 'module_perl'],
            ['combine_annovar_files', 'module_mugqic_tools']
        ],
        command="""\
/hpf/largeprojects/ccmbio/kng/mcgill_jacek/mugqic-2.2.0/utils/combine_annovar_files.pl \\
  --vcf {vcf_file} \\
  --exonvars {exonic_vars} \\
  --allvars {var_fn} \\
  --allvarsExtended {var_fn_extended} \\
  --table {table_file}{other_options} > {output}""".format(
            vcf_file=vcf_file,
            exonic_vars=exonic_vars,
            var_fn=var_fn,
            var_fn_extended=var_fn_extended,
            table_file=table_file,
            other_options=" \\\n  " + other_options if other_options else "",
            output=output
        ))


def convert2annovar(input, output, ini_section="convert2annovar"):
    
    other_options = config.param(ini_section, 'other_options', required=False)

    return Job(
        [input],
        [output],
        [
            ['convert2annovar', 'module_annovar'],
            ['convert2annovar', 'module_perl'],
        ],
        command="""\
convert2annovar.pl {input} \\
  --format vcf4 \\
  --includeinfo {other_options} > {output}""".format(
            input=input,
            other_options=other_options if other_options else "",
            output=output
    ))


def table_annovar(input, output, protocol, operation, ini_section="table_annovar"):
    
    annovar_db_dir = config.param(ini_section, 'annovar_db_dir', type='dirpath')
    other_options = config.param(ini_section, 'other_options', required=False)

    return Job(
        [input],
        [output + ".hg19_multianno.txt"],
        [
            ['table_annovar', 'module_annovar'],
            ['table_annovar', 'module_perl']
        ],
        command="""\
table_annovar.pl \\
  --buildver {buildver} \\
  --protocol {protocol} \\
  --operation {operation} \\
  --nastring . \\
  --outfile {output} \\
  --remove {other_options}\\
  {input} \\
  {annovar_db_dir}""".format(
            buildver=config.param(ini_section, 'assembly'),
            other_options=" \\\n  " + other_options if other_options else "",
            protocol=protocol,
            operation=operation,
            output=output,
            input=input,
            annovar_db_dir=annovar_db_dir
        ))
