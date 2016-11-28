#!/usr/bin/env python

# Python Standard Modules
import os

# MUGQIC Modules
from core.config import *
from core.job import *

def raw_counts(input_bam, output_dir):

    library_type = config.param('DEFAULT', 'library_type', type='string', required=True)  
    if library_type == 'single':
        strand = '--singleEnded'
    else:
        strand = '--stranded'
    

    return Job(
        [input_bam],
        [os.path.join(output_dir, 'QC.QORTS_COMPLETED_OK')],
        [['bedtools', 'module_bedtools']],
        command="""\
mkdir -p {output_dir} && \\
module load java/1.6.0 && \\
java -Xmx{ram} -jar {bfx}/QoRTs.jar QC {strand} {other_options} {input_bam} {gtf} {output_dir}""".format(
            input_bam = input_bam,
            output_dir = output_dir,
            strand = strand,
            other_options = config.param('jctseq_raw_counts', 'QoRTs_other_options', type='string', required=False),
            ram = config.param('jctseq_raw_counts', 'ram', type='string', required=True),
            gtf = config.param('jctseq_make_gff', 'gtf', type='filepath', required=True),
            bfx = config.param('DEFAULT', 'bfx_location', type='string', required=True),
            )
        )

def make_gff(output_gff):

    library_type = config.param('DEFAULT', 'library_type', type='string', required=True)
    if library_type == 'single':
        strand = ''
    else:
        strand = '--stranded'

    return Job(
        output_files=[output_gff],
        command="""\
mkdir -p jctseq && \\
module load java/1.6.0 && \\
java -Xmx{ram} -jar {bfx}/QoRTs.jar makeFlatGff {strand} {gtf} {output_gff}""".format(
            output_gff = output_gff,
            strand = strand,
            ram = config.param('jctseq_make_gff', 'ram', type='string', required=True),
            gtf = config.param('DEFAULT', 'gtf', type='filepath', required=True),
            bfx = config.param('DEFAULT', 'bfx_location', type='string', required=True)
            )
        )

def diff_prep(folder, contrast_name, gff):

    return Job(
        [gff],
        [os.path.join(folder, 'jscs1.r'), 'jctseq/jctseq.design'],
        [['bedtools', 'module_bedtools']],
        command="""\
cp {design_file} jctseq/jctseq.design ;\\
mkdir -p {folder} && \\
cp {bfx}/jscs1.r {folder} && \\
echo 'contrast <- decoder${contrast_name} \nfolder <- "{folder}"' >> {folder}/jscs1.r""".format(
            folder = folder,
            contrast_name = contrast_name,
            gff = gff,
            design_file = config.param('jctseq_diff_prep', 'design_file', type='filepath', required=True),
            bfx = config.param('DEFAULT', 'bfx_location', type='string', required=True)
            )
        )

def jscs_diff(input_gff, output_plots, output_results, jscs_file, raw_counts_list):

    return Job(
        raw_counts_list + [input_gff, jscs_file],
        [output_plots, output_results],
        [['bedtools', 'module_bedtools']],
        command="""\
mkdir -p {output_plots} && \\
mkdir -p {output_results} && \\
cat {bfx}/jscs2.r >> {jscs_file} && \\
module load R/3.3.0 && \\
module load gcc && \\
Rscript {jscs_file}""".format(
            input_gff = input_gff,
            output_plots = output_plots,
            output_results = output_results,
            jscs_file = jscs_file,
            raw_counts_list = raw_counts_list,
            bfx = config.param('DEFAULT', 'bfx_location', type='string', required=True)
            )
        )

def report(report_dependencies, report_file, report_template_dir, basename_report_file):

    event_names = config.param('miso_plot', 'events_names', type='string', required=True).split()

    events = ''
    for event in event_names:
        events += event + '.pdf\n'
    
    return Job(
        report_dependencies,
        [report_file],
        [['pandoc', 'module_pandoc']],
        command="""\
mkdir -p report && \\
pandoc --to=markdown \\
--template {report_template_dir}/{basename_report_file} \\
--variable events="{events}" \\
{report_template_dir}/{basename_report_file} \\
> {report_file}""".format(
            report_template_dir = report_template_dir,
            basename_report_file = basename_report_file,
            report_file = report_file,
            events = events
            ),
        report_files = [report_file],
        name = "jctseq_diff_report"
        )
