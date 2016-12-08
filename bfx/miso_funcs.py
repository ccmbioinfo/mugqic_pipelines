#!/usr/bin/env python

################################################################################

# New Miso step
#

################################################################################

# Python Standard Modules
import os

# MUGQIC Modules
from core.config import *
from core.job import *

def make_index(output_folder, gff, input_bam):

    var_set = '{print ".nochr"}'

    return Job(
        [input_bam],
        [output_folder],
        [['samtools', 'module_samtools']],
        command="""\
module load python/2.7.9 && \\
mkdir -p {output_folder} && \\
chr=$(samtools view  {input_bam} | head -n1 | awk '$3 !~ /^chr/ {var_set}') && \\
index_gff --index {gff}$chr.gff3 {output_folder}""".format(
            output_folder = output_folder,
            input_bam = input_bam,
            var_set = var_set,
            gff = gff
            )
        )

def paired_end(bam, sample_name):

    return Job(
        [bam],
        [os.path.join('miso', sample_name, 'insert-dist', sample_name + '.sorted.bam.insert_len')],
        [
            ['bedtools', 'module_bedtools'],
            ['samtools', 'module_samtools']
        ],
        command="""\
module load python/2.7.9 && \\
pe_utils --compute-insert-len {bam} {gff} --output-dir miso/{sample_name}/insert-dist/ --min-exon-size {min_exon_size} {other_options}""".format(
            bam = bam,
            sample_name = sample_name,
            min_exon_size = config.param('miso_paired_end', 'min_exon_size', type='string', required=True),
            gff = config.param('miso_paired_end', 'input_gff_file', type='filepath', required=True),
            other_options = config.param('miso_paired_end', 'other_options', type='string', required=False)
            )
        )

def compute_psi(input_bam, index_folder, output_header, output_folder_name, sample_name):
    
    p_end = ''
    insert_file = os.path.join('miso', sample_name,'insert-dist', sample_name + '.sorted.bam.insert_len')
    depend = []
    mean_com = 'BEGIN { FS=",";} /mean/ {print substr($1,7);}'
    sdev_com = 'BEGIN { FS=",";} /mean/ {print substr($2,6);}'
    len_com = '{print length($10)}'
 
    if config.param('miso_psi', 'library_type', type='string', required=True) == 'paired':
        depend = [insert_file]
        p_end = '--paired-end'

    return Job(
        [input_bam, index_folder] + depend,
        [output_header, output_folder_name],
        [
            ['python/2.7.9', 'module_python'],
            ['samtools', 'module_samtools']
        ],
        command="""\
module load python/2.7.9 && \\
samtools index {input_bam} && \\
mkdir -p {output_folder_name} && \\
unset mean sdev len && \\
mean=$(awk '{mean_com}' {insert_file}) ;\\
sdev=$(awk '{sdev_com}' {insert_file}) ;\\
len=$(samtools view {input_bam} | head -n1 | awk '{len_com}') && \\
miso --run {index_folder} {input_bam} --output-dir {output_folder_name} {other_options} --read-len $len {p_end} $mean $sdev""".format(
            input_bam=input_bam,
            index_folder=index_folder,
            output_header=output_header,
            output_folder_name=output_folder_name,
            p_end = p_end,
            mean_com = mean_com,
            sdev_com = sdev_com,
            len_com = len_com,
            insert_file = insert_file,
            other_options = config.param('miso_psi', 'other_options', type='string', required=False),
            )
        )

def summary(summarize_dir, summary):

    return Job(
        [summarize_dir],
        [summary],
        command="""\
module load python/2.7.9 && \\
summarize_miso --summarize-samples {summarize_dir} {summarize_dir} {other_options}""".format(
            summarize_dir=summarize_dir,
            summary=summary,
            other_options = config.param('miso_summarize', 'other_options', type='string', required=False)
            )
        )

def compare(sample_one, sample_two, output_directory, output, summaries_list):

    return Job(
        summaries_list,
        [output],
        command="""\
module load python/2.7.9 && \\
compare_miso --compare-samples miso/{sample_one} miso/{sample_two} {output_directory} {other_options}""".format(
            sample_one=sample_one,
            sample_two=sample_two,
            output=output,
            output_directory = output_directory,
            other_options = config.param('miso_diff', 'other_options', type='string', required=False)
            )
        )

def plot_settings(sample_names, bam_files, bfx):

    optional = ''
    logged = config.param('miso_plot_settings', 'logged', type='string', required=False)
    ymax = config.param('miso_plot_settings', 'ymax', type='string', required=False)
    sample_labels = config.param('miso_plot_settings', 'sample_labels', type='string', required=False)
    reverse_minus = config.param('miso_plot_settings', 'reverse_minus', type='string', required=False)
    nxticks = config.param('miso_plot_settings', 'nxticks', type='string', required=False)
    nyticks = config.param('miso_plot_settings', 'nyticks', type='string', required=False)
    bar_color = config.param('miso_plot_settings', 'bar_color', type='string', required=False)
    bf_thresholds = '' #config.param('miso_plot_settings', 'bf_thresholds', type='string', required=False),

    if logged != '':
        optional += 'logged = ' + logged + '\n'
    if ymax != '':
        optional += 'ymax = ' + ymax + '\n'
    if sample_labels != '':
        optional += 'sample_labels = ' + sample_labels + '\n'
    if reverse_minus != '':
        optional += 'reverse_minus = ' + reverse_minus + '\n'
    if nxticks != '':
        optional += 'nxticks = ' + nxticks + '\n'
    if nyticks != '':
        optional += 'nyticks = ' + nyticks + '\n'
    if bar_color != '':
        optional += 'bar_color = ' + bar_color + '\n'
    if bf_thresholds != '':
        optional += 'bf_thresholds = ' + bf_thresholds + '\n'

    return Job(
        output_files=['miso/sashimi_plot_settings.txt'],
        module_entries = [['samtools', 'module_samtools']],
        command="""\
mkdir -p miso/plots && \\
cp {bfx}/sashimi_plot_settings_temp.txt miso/sashimi_plot_settings.txt && \\
printf '\nbam_files = {bam_files} \nmiso_files = {sample_names} \n' >> miso/sashimi_plot_settings.txt && \\
printf '[plotting] \nfig_width = {fig_width} \nfig_height = {fig_height} \nintron_scale = {intron_scale} \nexon_scale = {exon_scale} \nshow_posteriors = {show_posteriors} \nbar_posteriors = {bar_posteriors} \ncolors = {colors} \ncoverages = {coverages} \n{optional}' >> miso/sashimi_plot_settings.txt""".format(
            bam_files = bam_files,
            sample_names = sample_names,
            bfx = bfx,
            fig_width = config.param('miso_plot_settings', 'fig_width', type='string', required=True),
            fig_height = config.param('miso_plot_settings', 'fig_height', type='string', required=True),
            intron_scale = config.param('miso_plot_settings', 'intron_scale', type='string', required=True),
            exon_scale = config.param('miso_plot_settings', 'exon_scale', type='string', required=True),
            show_posteriors = config.param('miso_plot_settings', 'show_posteriors', type='string', required=True),
            bar_posteriors = config.param('miso_plot_settings', 'bar_posteriors', type='string', required=True),
            colors = config.param('miso_plot_settings', 'colors', type='string', required=True),
            coverages = config.param('miso_plot_settings', 'coverages', type='string', required=True),
            optional = optional
            )
        )

def plot_insert_len(insert_file):

    return Job(
        [insert_file, 'miso/sashimi_plot_settings.txt'],
        ['miso/plots'],
        [['samtools', 'module_samtools']],
        command="""\
module load python/2.7.9 && \\
sashimi_plot --plot-insert-len {insert_file} miso/sashimi_plot_settings.txt  --output-dir miso/plots""".format(
            insert_file = insert_file
            )
        )

def plot_bf_dist(bf_file):

    return Job(
        [bf_file, 'miso/sashimi_plot_settings.txt'],
        ['miso/plots'],
        [['samtools', 'module_samtools']],
        command="""\
module load python/2.7.9 && \\
sashimi_plot --plot-bf-dist {bf_file} miso/sashimi_plot_settings.txt --output-dir miso/plots""".format(
            bf_file = bf_file
            )
        )

def plot_event(dependencies, event_name):
    
    output_plot = os.path.join('miso', 'plots', event_name + '.pdf')

    return Job(
        dependencies + ['miso/indexed_SE', 'miso/sashimi_plot_settings.txt'],
        [output_plot],
        [['samtools', 'module_samtools']],
        command="""\
module load python/2.7.9 && \\
sashimi_plot --plot-event "{event_name}" miso/indexed miso/sashimi_plot_settings.txt --output-dir miso/plots""".format(
            event_name = event_name,
            output_plot = output_plot
            )
        )

def report(dependencies, report_file, report_template_dir, basename_report_file):

    plots = ''
    events = ''
    events_names = config.param('miso_plot', 'events_names', type='string', required=True).split()
    for event_name in events_names:
        event_plot = os.path.join('miso', 'plots', event_name + '.pdf')
        
        dependencies.append(event_plot)

        plots += event_plot + ' '

        events += event_name + '.pdf\n'

    return Job(
        dependencies,
        [report_file],
        [['pandoc', 'module_pandoc']],
        command="""\
mkdir -p report && \\
cp {plots}report && \\
pandoc --to=markdown \\
--template {report_template_dir}/{basename_report_file} \\
--variable events="{events}" \\
{report_template_dir}/{basename_report_file} \\
> {report_file}""".format(
            report_template_dir = report_template_dir,
            basename_report_file = basename_report_file,
            report_file = report_file,
            events = events,
            plots = plots,
            ),
        report_files = [report_file],
        name = "miso_plot_report"
        )
