#!/usr/bin/env python

# Python Standard Modules
import os

# MUGQIC Modules
from core.config import *
from core.job import *

def make_fastq(input_bam, output_fastq):
    return Job(
        [input_bam],
        [output_fastq],
        [['bedtools', 'module_bedtools']],
        command="""\
module load tophat && \\
bam2fastx --fastq -M -o {output_fastq} -N {input_bam}""".format(
            input_bam = input_bam,
            output_fastq = output_fastq,
            )
        )


def alignment(input_fastq, align_output):
    return Job(
        [input_fastq],
        [align_output],
        [['bedtools', 'module_bedtools']],
        command="""\
module load vast-tools && \\
mkdir -p vast_out && \\
vast-tools align {input_fastq} --output vast_out {other_options} -sp {species}""".format(
            input_fastq = input_fastq,
            align_output = align_output,
            other_options = config.param('vast_tools_align', 'other_options', type='string', required=False),
            species = config.param('vast_tools_align', 'species', type='string', required=True)
            )
        )

def comb(dependency_list, inclusion_table):
    return Job(
        dependency_list,
        [inclusion_table],
        [['bedtools', 'module_bedtools']],
        command="""\
module load vast-tools && \\
vast-tools combine -o vast_out {other_options} -sp {species}""".format(
            dependency_list = dependency_list,
            inclusion_table = inclusion_table,
            other_options = config.param('vast_tools_combine', 'other_options', type='string', required=False),
            species = config.param('vast_tools_align', 'species', type='string', required=True)
            )
        )


def differential_splicing(sample_one, sample_two, input_table, inclusion_table):
    return Job(
        [input_table],
        [inclusion_table],
        [['bedtools', 'module_bedtools']],
        command="""\
module load vast-tools && \\
vast-tools diff -a {sample_one} -b {sample_two} -o vast_out {other_options} > {inclusion_table}""".format(
            sample_one = sample_one,
            sample_two = sample_two,
            input_table = input_table,
            other_options = config.param('vast_tools_diff', 'other_options', type='string', required=False),
            inclusion_table = inclusion_table,            
            )
        )

def plots(inclusion_table, events_list):

    search_str = '/EVENT/ || '
    for event in events_list:
        search_str += '/' + event + '/ || '

    return Job(
        [inclusion_table],
        ['vast_out/plot_events.tab', 'vast_out/plot_events.PSI_plots.pdf'],
        [['bedtools', 'module_bedtools']],
        command="""\
module load vast-tools && \\
awk '{search}' {inclusion_table} > vast_out/plot_events.tab && \\
vast-tools plot {other_options} vast_out/plot_events.tab""".format(
            inclusion_table = inclusion_table,
            search = search_str[0:-4],
            other_options = config.param('vast_tools_plot', 'other_options', type='string', required=False)
            )
        )

def report(report_file, report_template_dir, basename_report_file, full_inclusion):

    return Job(
        [full_inclusion, 'vast_out/plot_events.PSI_plots.pdf'],
        [report_file],
        [['pandoc', 'module_pandoc']],
        command="""\
mkdir -p report && \\
cp {full_inclusion} report/INCLUSION_LEVELS_FULL.tab && \\
cp vast_out/plot_events.PSI_plots.pdf report && \\
pandoc --to=markdown \\
--template {report_template_dir}/{basename_report_file} \\
{report_template_dir}/{basename_report_file} \\
> {report_file}""".format(
            report_template_dir = report_template_dir,
            basename_report_file = basename_report_file,
            report_file = report_file,
            full_inclusion = full_inclusion
            ),
        report_files = [report_file],
        name = "vast_tools_report"
        )

