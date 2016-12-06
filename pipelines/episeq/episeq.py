#!/usr/bin/env python
"""
Epigenetics pipeline for RRBS/WGBS data
=======================================

Background
----------
EpiSeq is a differential analysis pipeline for BS-seq sequencing. Currently only RRBS and WGBS datasets are
tested to work with this pipeline. Similar to the other MUGQIC pipeline series, EpiSeq uses two metadata files to
set up the pipeline. The design file is used to group samples into case vs control. The readsets files refer to
which input files correspond to each same and, if applicable, which pairs of input files correspond to paired reads.
A readset is considered to be one instance/lane/segment of sequencing data. This is often used when multiple libraries
are used for a given sample, or if multiplexing was done. These techniques tend to generate multiple sets of data
for a given sample. The readset file allows users to specify these relationships.

This simple pipeline uses Bismark, Trim Galore!, R, and Picard to process data. See the README file for more information
about the steps in the pipeline.

Input
-----
- ``FASTQ`` or ``BAM`` files containing methylation sequencing data
- Reference genome in ``FASTA`` format
- MUGQIC formatted ``.design`` file
- MUGQIC formatted ``.readset`` file
- Episeq pipeline's ``.ini`` file
- Episeq pipeline's ``bismark_merge_reports.py`` script

Output Data
-----------
- Bam file:
    - dedup/{sample.name}/{sample.name}.merged.deduplicated.bam
- Methylation Graph Plot:
    - methyl_calls,{sample.name}/{sample.name}.merged.deduplicated.bedGraph.gz
- Methylation Coverage:
    - methyl_calls{sample.name}/{sample.name}.merged.deduplicated.bismark.cov.gz
- CpG Methylation List:
    - methyl_calls/{sample.name}/{sample.name}.merged.deduplicated.CpG_report.txt.gz
- Differentially Methylated Positions:
    - differential_methylated_positions/{contrast.name}_RRBS_differential_methylated_pos.csv
- Differentially Methylated Regions:
    - differential_methylated_regions/{contrast.name}_RRBS_differential_methylated_regions.csv

Output Reports
--------------
- Alignment report:
    - merged/{sample.name}/{sample.name}.merged_aligned_PE_report.txt
- Deduplication report:
    - dedup/{sample.name}/{sample.name}.merged.deduplication_report.txt
- Nucleotide Coverage:
    - dedup/{sample.name}/{sample.name}.merged.deduplicated.nucleotide_stats.txt
    - merged/{sample.name}/{sample.name}.merged.nucleotide_stats.txt
- Methylation Extraction Report:
    - methyl_calls/{sample.name}/{sample.name}.merged.deduplicated_splitting_report.txt
- M-bias Plot:
    - methyl_calls/{sample.name}/{sample.name}.merged.deduplicated.M-bias.txt
- Bismark HTML Sample Summary:
    - bismark_summary_report/{sample.name}_final_bismark_report.html

Workflow
--------
Below is the steps defined in this pipeline. The alphabetical designation shows which steps run concurrently in the
pipeline. It shows the order of the pipeline and the dependencies for each step.

 1. (a) Genome Methylation Conversion (bismark_prepare_genome)
 2. (a) Pre-Trim Quality Check (pre_qc_check)
 3. (a) Read and Adapter Trimming (trim_galore)
 4. (b) Alignment (bismark_align)
 5. (c) Merge Alignment Statistics(merge_bismark_alignment_report)
 6. (c) Merge Aliged BAMs (picard_merge_sam_files)
 7. (d) Recalculate nucleotide coverage (merged_nuc_stats)
 8. (d) Remove Duplicate Reads (bismark_deduplication)
 9. (e) Recalculate nucleotide coverage 2 (calc_dedup_nucleotide_coverage)
10. (e) Methylation Calling and Analysis (bismark_methylation_caller)
11. (f) Bismark Sample-Level Report Generator (bismark_html_report_generator)
12. (g) Position Specific Differential Analysis (differential_methylated_pos)
13. (g) Regional Differential Analysis (differential_methylated_regions)

Implementation
--------------
The pipeline is arranged to gradually decrease the number of input files when it make sense to do so. The trimming
and QC step would be too slow if we merge readsets together, so we decided to trim each file individually. Then,
the alignment phase allows us to convert out fastQ files to BAM files. This allows us to merge paired datasets into
one BAM file per readset without explicitly doing so. After this, readsets are combined to run the analysis steps.
"""

# Python Standard Modules
import argparse
import collections
import logging
import os
import re
import sys

# Append mugqic_pipelines directory to Python library path
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0])))))

# MUGQIC Modules
from core.config import *
from core.job import *
from core.pipeline import *
from bfx.readset import *
from bfx.design import *

from bfx import rmarkdown
from pipelines import common
import utils

log = logging.getLogger(__name__)


class Episeq(common.Illumina):
    """
    The Episeq pipeline takes FASTQ or BAM files (unsorted) as input as well as two metadata files and a configuration
    file. Refer to the user guide for more information on running the pipeline.
    """

    def __init__(self):
        self.argparser.add_argument("-d", "--design", help="design file", type=file)
        super(Episeq, self).__init__()

    @staticmethod
    def bismark_prepare_genome():
        """
        Bismark requires a processed reference genome to compare with the epigenome. This step can take several hours,
        depending on the size of the genome. It creates an index of bisulfite conversions and takes make space than
        the genome itself. This module will only create the output in the same directory as the genome file, so
        a symlink is needed to create a "copy" in the desired output directory. Only runs once and generates one job.

        Note: The genome should be in fasta format.

        Input: A reference sequence file as specified by user. Configuration is set in the episeq.ini file.
        Output: A directory called Bisulfite_Genome and a dinucleotide composition report.

        :return: A list of jobs that needs to be executed in this step.
        :rtype: list(Job)
        """
        # Variable preparations
        ref_seq = config.param('bismark_prepare_genome', 'genome_file', type='filepath')
        local_ref_seq = os.path.join('bismark_prepare_genome', os.path.basename(ref_seq))
        output_idx = "bismark_prepare_genome/Bisulfite_Genome"
        modules = [['bismark_prepare_genome', 'module_bowtie2'],
                   ['bismark_prepare_genome', 'module_samtools'],
                   ['bismark_prepare_genome', 'module_perl'],
                   ['bismark_prepare_genome', 'module_bismark']]

        # Job creation
        mkdir_job = Job(command='mkdir -p bismark_prepare_genome')
        link_job = Job(input_files=[ref_seq],
                       command='cp -sfu ' + os.path.abspath(ref_seq) + ' ' + os.path.abspath(local_ref_seq),
                       removable_files=[local_ref_seq])
        main_job = Job(output_files=[output_idx], module_entries=modules,
                       command="bismark_genome_preparation --verbose bismark_prepare_genome/")
        nuc_count = Job(output_files=[], module_entries=modules,
                        command="bam2nuc --genomic_composition_only --genome_dir bismark_prepare_genome/ "
                                "--dir bismark_prepare_genome/")
        return [concat_jobs([mkdir_job, link_job, main_job, nuc_count],
                            name='bismark_prepare_genome.' + os.path.basename(ref_seq))]

    def pre_qc_check(self):
        """
        Runs FastQC on the unprocessed fastq files to generate a baseline report. Helpful when comparing to post-trim
        metrics.

        Note: While FastQC does use perl, do not load any perl modules until the sha-bang statement uses /usr/bin/env.
        Failing to do so will result in incorrect libraries being used as /usr/bin/perl is not likely the same version
        as what you load.

        A single report file is generated with the help of a intermediate file
        that lets tasks update the table format as they complete. The drawback
        is that pandoc will have to run for each task, but is to allow some
        jobs to fail and still have a working report.

        :return: A list of jobs that needs to be executed in this step.
        :rtype: list(Job)
        """
        # static variables
        jobs = []  # List of jobs to run
        requires = [['pre_qc_check', 'module_java'],
                    ['pre_qc_check', 'module_fastqc']]
        tmpdir = config.param('pre_qc_check', 'temp_dir', required=False) or config.param('DEFAULT', 'tmp_dir')
        template_var_hold = 'pre_qc_check/' + datetime.datetime.now().strftime('%Y_%m_%d') + "_template_var_strings.txt"
        report_file = 'report/EpiSeq.pre_qc_check.md'
        report_data = 'report/data/pre_qc_check'
        fill_in_templ = "| {sample} | {readset} | [Report 1]({fq1_report}) [Download 1]({fq1_download}) | " \
                        "[Report 2]({fq2_report}) [Download 2]({fq2_download}) | {submission} | {completion} |"
        # Generate processing job
        for sample in self.samples:
            for num, readset in enumerate(sample.readsets, 1):
                # Determine what fastq files are given. 0 <= len(raw_fq) <= 2
                raw_fq = filter(None, [readset.fastq1, readset.fastq2])
                if not raw_fq:
                    # If no fastq files, the sample just might only have a bam file, else throw error
                    log.info('No fastq files found for readset for ' + readset.name + '. Skipping...')
                    if readset.bam:
                        continue
                    else:
                        raise ValueError('No sequencing data found for ' + readset.name)
                # Determine output directory format
                # My method of collapsing unneeded directories. If we only have 1 readset, don't make nested directories
                if len(sample.readsets) == 1:
                    out_dir = os.path.join('pre_qc_check', sample.name)
                else:
                    out_dir = os.path.join('pre_qc_check', sample.name, readset.name)
                # Generate output path
                # Obtain the output of the filename. ex. 'SRS23542_1' + '...'
                id_name = [os.path.basename(nom).split('.gz')[0].split('.fastq')[0] for nom in raw_fq]
                file_names = [nom + ext for nom in id_name for ext in ['.html', '.zip']]
                output = [os.path.join(out_dir, name) for name in file_names]
                # Job creation
                mkdir_job = Job(command='mkdir -p ' + out_dir)
                job = Job(input_files=raw_fq,
                          output_files=output,
                          module_entries=requires,
                          command="fastqc -o {out_dir} {others} -d {tmpdir} {inputs}".format(
                              inputs=' '.join(raw_fq),
                              out_dir=out_dir,
                              others=config.param('pre_qc_check', 'other_options', required=False),
                              tmpdir=tmpdir),
                          removable_files=[output])
                # Generate new report with new entry, only way to ensure nothing is missed, given no order guarantee.
                # Since order is mixed, we can't make an organized table for viewing (without some fancier code)
                report_body = fill_in_templ.format(
                    sample=sample.name,
                    readset=readset.name,
                    fq1_report=os.path.join(report_data, id_name[0]),
                    fq1_download=os.path.join(report_data, id_name[0] + '.zip'),
                    fq2_report=os.path.join(report_data, id_name[1]) if len(id_name) == 2 else '',
                    fq2_download=os.path.join(report_data, id_name[1] + '.zip') if len(id_name) == 2 else '',
                    submission=datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
                    completion="$(date +%Y-%m-%d %H:%M:%S)"
                )
                report_body += '\n'  # End entry

                # Code to fill in template.
                command = """\
# Mutual exclusion with other jobs
flock -x ".{table_hold}.lock" -c "echo -e {entry} >> {table_hold}" && \\
mkdir -p {data_loc} && \\
cp -f {output_file} {data_loc} && \\
table=$(cat {table_hold}) && \\
pandoc \\
  {report_template_dir}/{basename_report_file} \\
  --template {report_template_dir}/{basename_report_file} \\
  --variable data_table="$table" \\
  --to markdown \\
  > {report_file}""".format(
                    entry=report_body,
                    report_template_dir=self.report_template_dir,
                    basename_report_file=os.path.basename(report_file),
                    output_file=" ".join([os.path.join(report_data, os.path.basename(out_log)) for out_log in output]),
                    data_loc=report_data,
                    table_hold=template_var_hold,
                    report_file=report_file)
                update_template = Job(
                    output_files=[os.path.join(report_data, os.path.basename(out_log)) for out_log in output],
                    module_entries=[['pre_qc_check', 'module_pandoc']],
                    command=command,
                    report_files=[report_file])
                # Add to list of jobs
                jobs.append(concat_jobs([mkdir_job, job, update_template],
                                        name='pre_qc_check.' + readset.name))
        return jobs

    def trim_galore(self):
        """
        This step trims raw FASTQ files for quality control using Trim Galore!
        This is a pre-processing step to ensure quality control.

        To run Trim Galore in paired mode, two fastq files must be specified in the readset file with the following
        pairewise naming convention: file1_1.fq file1_2.fq SRR2_1.fq.gz SRR2_2.fq.gz

        Note: Do not load perl modules because Trim Galore and FastQC only uses /usr/bin/perl. It will not use the
        correct perl binary and thus the libraries are set up for the wrong version.

        :return: A list of jobs that needs to be executed in this step.
        :rtype: list(Job)
        """
        # Report Prep
        report_file = 'report/EpiSeq.trim_galore.md'
        report_data = 'report/data/trim_galore'
        template_string_file = 'trimmed/' + datetime.datetime.now().strftime('%Y_%m_%d') + "_template_var_strings.txt"
        fill_in_entry = '| {sample} | {readset} | {trim1_view} {trim1_download} {trim2_view} {trim2_download} | ' \
                        '{qc1_view} {qc1_download} {qc2_view} {qc2_download} | {submission} | {completion} |'
        # Awk script to parse trim reports view easy viewing
        # It will read the start of each line to determine if it is something that we want to keep, or discard.
        # Headers have to be reformatted to a lower level, making the document make sense within the template.
        # Additional code is added to properly add table syntax when extra blank lines are found.
        output_parser = """\
BEGIN {table=0;}
{
    if ($0 ~ /^SUMMARISING/) {
        print "##### " $0;
        table=0;
    } else if ($0 ~ /^RUN STATISTICS/) {
        print "###### " $0;
        table=0;
    } else if ($0 ~ /^=+$/) {
        table=0;
        next;
    } else if (match($0, /^=== ([[:alnum:][:blank:]]+) ===$/, a)) {
        print "###### " a[1];
        table=0;
    } else if ($0 ~ /^\|Input filename.+/) {
        print ""; print "|Attribute|Value|";
        print "|----|----|"; print $0;
        table=1;
    } else if ($0 ~ /^\|Total reads.+/) {
        print "";
        print "|Metric|Value|";
        print "|----|----|"; print $0;
        table=1;
    } else if ($0 ~ /^length\|.+/) {
        print "";
        print $0;
        print "----:|----:|----:|----:|----:";
        table=1;
    } else if ($0 ~ /^Bases preceding.+/) {
        print $0;
        print "";
        print "|Base|Value|";
        print "|:----:|----:|";
        table=1;
    } else if ($0 ~ /^Running/ || $0 ~/^Output/ || $0 ~ /^This is cut.+/) {
        next;
    } else if (/^\|/) {
        print $0;
        table=1;
    } else if ($0 ~ /^$/ && table == 1) {
        next;
    } else if ($0 ~ /^$/ && table == 0) {
        print $0;
    } else {
        print $0;
    }
}"""
        jobs = []
        for sample in self.samples:
            for readset in sample.readsets:  # Iterate through each readset in project
                # Determine if we can run trimming step for this set of data (in case user skips all other steps)
                if readset.bam and not filter(None, [readset.fastq1, readset.fastq2]):
                    log.info('No fastq files found for readset for ' + readset.name + '. Skipping...')
                    continue  # There's still a bam file to process later, so skip this for now
                elif not (readset.bam or readset.fastq1 or readset.fastq2):
                    raise ValueError("There are no data associated with this readset: " + readset.name + "!\n")

                if len(sample.readsets) == 1:  # output directory, colapse when needed
                    trim_directory = os.path.join("trimmed", sample.name)
                else:
                    trim_directory = os.path.join("trimmed", sample.name, readset.name)
                run_type = readset.run_type
                protocol = readset.library
                input_files = filter(None, [readset.fastq1, readset.fastq2])
                file_basename = [os.path.join(trim_directory, os.path.basename(in_file).split('.')[0])
                                 for in_file in input_files]

                # List all possible outputs
                if len(input_files) == 2:
                    input1_logs = [trim_directory + '/' + os.path.basename(readset.fastq1) + '_trimming_report.txt',
                                   file_basename[0] + "_val_1_fastqc.html", file_basename[0] + "_val_1_fastqc.zip",
                                   file_basename[0] + "_val_1.fq.gz"]
                    input2_logs = [trim_directory + '/' + os.path.basename(readset.fastq2) + '_trimming_report.txt',
                                   file_basename[1] + "_val_2_fastqc.html", file_basename[1] + "_val_2_fastqc.zip",
                                   file_basename[1] + "_val_2.fq.gz"]
                else:
                    input1_logs = [trim_directory + '/' + os.path.basename(readset.fastq1) + '_trimming_report.txt',
                                   file_basename[0] + "_trimmed_fastqc.html", file_basename[0] + "_trimmed_fastqc.zip",
                                   file_basename[0] + "_trimmed.fq.gz"]
                    input2_logs = []

                # Parse custom options that may affect output (ex. Running FastQC or suppressing reports)
                add_options = config.param('trim_galore', 'other_options').split()
                run_qc = not ('--fastqc_args' in add_options or '--fastqc' in add_options)
                report_out = '--no_report_file' not in add_options

                # Trim Galore has no built in option to change the filenames of the output
                # Below are the default output names when running in paired or single mode
                logs = list()
                if run_type == "PAIRED_END":
                    if not readset.fastq2:
                        raise ValueError("Expecting paired reads be named as follows file1_1.fq file1_2.fq.")
                    if report_out:
                        logs += [input1_logs[0]] + [input2_logs[0]]  # Trim report
                    if run_qc:
                        logs += [input1_logs[1]] + [input2_logs[1]]  # FastQC html
                        logs += [input1_logs[2]] + [input2_logs[2]]  # FastQC zip
                        output_files = [input1_logs[3]] + [input2_logs[3]]
                else:
                    if report_out:
                        logs += [input1_logs[0]]
                    if run_qc:
                        logs += [input1_logs[1]] + [input1_logs[2]]
                    output_files = [input1_logs[3]]

                # Define jobs
                mkdir_job = Job(command="mkdir -p " + trim_directory)
                trim_job = Job(input_files=input_files,
                               output_files=logs + output_files,
                               module_entries=[['trim_galore', 'module_fastqc'],
                                               ['trim_galore', 'module_java'],
                                               ['trim_galore', 'module_trim_galore'],
                                               ['trim_galore', 'module_cutadapt']],
                               command="trim_galore {protocol} {library_type} {other} {directory} {fastq}".format(
                                   library_type="--paired" if run_type == "PAIRED_END" else "",
                                   protocol='--rrbs' if protocol == 'RRBS' else '',
                                   other=config.param("trim_galore", "other_options"),
                                   directory='--output_dir ' + trim_directory,
                                   fastq=' '.join(input_files)),
                               removable_files=output_files)
                entry = fill_in_entry.format(sample=sample.name,
                                             readset=readset.name,
                                             trim1_view='[View details 1](' +
                                                        os.path.join(report_data,
                                                                     os.path.basename(input1_logs[0])) + '.md)',
                                             trim1_download='[Download details 1](' +
                                                            os.path.join(report_data,
                                                                         os.path.basename(input1_logs[0])) + ')',
                                             trim2_view='[View details 2](' +
                                                        os.path.join(report_data,
                                                                     os.path.basename(input2_logs[0])) + '.md)'
                                             if input2_logs else '',
                                             trim2_download='[Download details 2](' +
                                                            os.path.join(report_data,
                                                                         os.path.basename(input2_logs[0])) + ')'
                                             if input2_logs else '',
                                             qc1_view='[View report 1](' +
                                                      os.path.join(report_data,
                                                                   os.path.basename(input1_logs[1])) + ')',
                                             qc1_download='[Download report 1](' +
                                                          os.path.join(report_data,
                                                                       os.path.basename(input1_logs[2])) + ')',
                                             qc2_view='[View report 2](' +
                                                      os.path.join(report_data,
                                                                   os.path.basename(input2_logs[1])) + ')'
                                             if input2_logs else '',
                                             qc2_download='[Download report 2](' +
                                                          os.path.join(report_data,
                                                                       os.path.basename(input2_logs[2])) + ')'
                                             if input2_logs else '',
                                             submission=datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
                                             completion="$(date +%Y-%m-%d %H:%M:%S)")
                entry += '\n'  # end entry
                new_logs = [os.path.join(report_data, os.path.basename(out_log)) for out_log in
                            filter(None, input1_logs + input2_logs)]
                command = """\
    # Mutual exclusion with other jobs
    flock -x ".{table_hold}.lock" -c "echo -e {entry} >> {table_hold}" && \\
    mkdir -p {data_loc} && \\
    cp -f {output_file} {data_loc} && \\
    for i in {individual_page}; do
        sed -r 's%^([^0-9S][a-z ].+):(\s+.+)%|\\1|\\2|%g' $i | sed 's/\t/|/g' | awk '{script}' > "$i.md"
    done
    table=$(cat {table_hold}) && \\
    pandoc \\
    {report_template_dir}/{basename_report_file} \\
        --template {report_template_dir}/{basename_report_file} \\
        --variable data_table="$table" \\
        --to markdown \\
    > {report_file}""".format(
                    entry=entry,
                    report_template_dir=self.report_template_dir,
                    basename_report_file=os.path.basename(report_file),
                    output_file=" ".join(logs),
                    individual_page=new_logs,
                    data_loc=report_data,
                    table_hold=template_string_file,
                    script=output_parser,
                    report_file=report_file)
                report_job = Job(
                    output_files=[report_file] + new_logs,
                    command=command,
                    module_entries=[['trim_galore', 'module_pandoc']],
                    report_files=[report_file])
                jobs.append(concat_jobs([mkdir_job, trim_job, report_job], name='trim_galore.' + readset.name))
            return jobs

    def bismark_align(self):
        """
        This step aligns trimmed reads to a bisulfite converted reference genome using Bismark. This create
        BAM files and will only be compressed if the input is also compressed (which usually will be the case).
        All readsets are aligned individually, but combines paired files into one output file. The output files
        are all placed in the same sample directory. No sub-directories are made.

        This step requires bismark_prepare_genome and the relevant trim_galore step.

        Note: Despite what the manual says, the source code shows that -un and --ambiguous produces fq files, not txt.

        Input: Trimmed version of input files as a fastq file. (trimmed/*)
        Output: A BAM/SAM file in aligned/<sample_name>/*

        :return: A list of jobs that needs to be executed in this step.
        :rtype: list(Job)
        """
        # Report Generation Setup
        report_file = 'report/EpiSeq.bismark_align.md'
        report_data = 'report/data/bismark_align'
        template_string_file = 'aligned/' + datetime.datetime.now().strftime('%Y_%m_%d') + "_template_var_strings.txt"
        fill_in_entry = '| {sample} | {readset} | {report_view} | {coverage_view} | {submission} | {completion} |'
        awk_script = """\
BEGIN { table=0; };
{
    if ($0 ~ /^$/) {
        if (table==1) {
            next;}
        else {
            print $0;
        }
    }
    else if ($0 ~ /^[^|]+|[^|]+$/ && table==0) {
            print ""; print "|Item|Value|"; print "|-----|----:|"; print $0;
            table=1;
    }
    else if($0 ~ /^Final.+/) {
            table=0;
            print "";
            print "##### " $0;
    }
    else if($0 ~ /^Number of sequence.+:$/) {
        table=0''
        print "";
        print $0;
    }
    else {
        print $0;
    }
}"""
        jobs = []
        for sample in self.samples:  # Process samples
            for readset in sample.readsets:
                # Special cases to consider when finding files from the trimming step
                if len(sample.readsets) == 1:
                    trim_prefix = os.path.join("trimmed", sample.name)
                else:
                    trim_prefix = os.path.join("trimmed", sample.name, readset.name)

                # Common Variables
                run_type = readset.run_type
                align_directory = os.path.join("aligned", sample.name)
                output_basename = os.path.join(align_directory, readset.name)
                user_options = config.param('bismark_align', 'other_options').split()
                input_basename = [os.path.join(trim_prefix, os.path.basename(in_file).split('.')[0])
                                  for in_file in filter(None, [readset.fastq1, readset.fastq2])]

                # Check what input files are found
                if not input_basename:
                    if readset.bam:
                        continue
                    else:
                        raise IOError("""{readset} has no input files!""".format(readset=readset.name))

                # Again, the suffix is hardcoded into the script. So we have to match it too. PE and SE have diff names
                if run_type == "PAIRED_END" and readset.fastq2:
                    input_files = [input_basename[0] + '_val_1.fq.gz', input_basename[1] + "_val_2.fq.gz"]
                    cmd_in = '-1 {fastq1} -2 {fastq2}'.format(fastq1=input_files[0], fastq2=input_files[1])
                    out_files = [output_basename + "_aligned_pe.bam"]
                    report_log = [output_basename + "_aligned_PE_report.txt"]
                    # Optional output files, depending on flags specified
                    if '--nucleotide_coverage' in user_options:
                        report_log += [output_basename + "_aligned_pe.nucleotide_stats.txt"]
                    if '-un' in user_options or '--unmapped' in user_options:
                        out_files += [output_basename + "_aligned_unmapped_reads_1.fq.gz",
                                      output_basename + "_aligned_unmapped_reads_2.fq.gz"]
                    if '--ambiguous' in user_options:
                        out_files += [output_basename + "_aligned_ambiguous_reads_1.fq.gz",
                                      output_basename + "_aligned_ambiguous_reads_2.fq.gz"]
                    if '--ambig_bam' in user_options:
                        out_files += [output_basename + '_aligned_pe.ambig.bam']
                elif run_type == "SINGLE_END":
                    input_files = [input_basename[0] + "_trimmed.fq.gz"]
                    cmd_in = '--single_end {fastq1}'.format(fastq1=input_files[0])
                    out_files = [output_basename + "_aligned.bam"]
                    report_log = [output_basename + "_aligned_SE_report.txt"]
                    # Optional output files - depends on flag specified.
                    if '--nucleotide_coverage' in user_options:
                        report_log += [output_basename + "_aligned.nucleotide_stats.txt"]
                    if '-un' in user_options or '--unmapped' in user_options:
                        out_files += [output_basename + "_unmapped_reads.fq.gz"]
                    if '--ambiguous' in user_options:
                        out_files += [output_basename + "_ambiguous_reads.fq.gz"]
                    if '--ambig_bam' in user_options:
                        out_files += [output_basename + '_aligned.ambig.bam']
                else:
                    raise AttributeError("Unknown run_type: " + run_type + ". Unknown file output name.")

                # Job creation
                mkdir_job = Job(command="mkdir -p " + align_directory)
                job = Job(
                    input_files + ["bismark_prepare_genome/Bisulfite_Genome"],
                    out_files + report_log,
                    [["bismark_align", "module_bowtie2"],
                     ["bismark_align", "module_samtools"],
                     ['bismark_align', 'module_perl'],
                     ['bismark_align', 'module_bismark']],
                    command="""\
bismark -q {other} --temp_dir {tmpdir} --output_dir {directory} \
    --basename {basename} --genome_folder bismark_prepare_genome {input}
        """.format(
                        directory=align_directory,
                        other=config.param("bismark_align", "other_options"),
                        tmpdir=config.param('bismark_align', 'tmp_dir', required=False) or
                            config.param('DEFAULT', 'tmp_dir', required='True'),
                        input=cmd_in,
                        basename=readset.name + '_aligned'),
                    removable_files=out_files
                )

                # Generate report stub
                new_logs = [os.path.join(report_data, os.path.basename(txt)) for txt in report_log]
                report_entry = fill_in_entry.format(
                    sample=sample.name,
                    readset=readset.name,
                    report_view='[View HTML]({a}) [Download Raw]({b})'.format(
                        a=os.path.join(report_data, report_log[0] + '.md'),
                        b=os.path.join(report_data, report_log[0])),
                    coverage_view='[View HTML]({a}) [Download Text]({b})'.format(
                        a=os.path.join(report_data, report_log[1] + '.md'),
                        b=os.path.join(report_data, report_log[1])) if report_log[1] else "N/A",
                    submission=datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
                    completion="$(date +%Y-%m-%d %H:%M:%S)")
                command = """\
# Mutual exclusion with other jobs
flock -x .{table_hold}.lock -c "echo -e {entry} >> {table_hold}" && \\
mkdir -p {data_dir} && \\
cp -f {reports} {data_dir} && \
for i in {logs}; do
    sed -r 's/:\s+/|/g' $i | egrep -v "^Option" | egrep -v "^=+$" | awk '{script}' > $i".md";
done
table=$(cat {table_hold}) && \\
pandoc \\
{report_template_dir}/{basename_report_file} \\
    --template {report_template_dir}/{basename_report_file} \\
    --variable data_table="$table" \\
    --to markdown \\
> {report_file}
""".format(
                    table_hold=template_string_file,
                    entry=report_entry,
                    data_dir=report_data,
                    reports=' '.join(report_log),
                    logs=' '.join(new_logs),
                    script=awk_script,
                    report_template_dir=self.report_template_dir,
                    basename_report_file=os.path.basename(report_file),
                    report_file=report_file
                )
                report_job = Job(output_files=[report_file] + new_logs,
                                 report_files=[report_file],
                                 module_entries=[['bismark_align', 'module_pandoc']],
                                 command=command)
                jobs.append(concat_jobs([mkdir_job, job, report_job], name="bismark_align." + readset.name))
        return jobs

    def merge_bismark_alignment_report(self):
        """
        This steps takes all of Bismark's alignment reports for a sample and merges them with a custom script.
        Some stats are recalculated to match the total read population and specific settings are lost due to the
        aggregation. Outputs the file to the merge directory, which will contain other merged results and outputs.

        :return: A list of jobs that needs to be executed in this step.
        :rtype: list(Job)
        """
        jobs = []
        for sample in self.samples:
            log_reports = []
            align_directory = os.path.join("aligned", sample.name)  # Previous step's output dir

            # Get all log reports for this sample
            for readset in sample.readsets:
                if readset.bam or not (readset.fastq1 or readset.fastq2):
                    continue
                # Get correct name values
                log_basename = os.path.join(align_directory, readset.name)
                if readset.run_type == "PAIRED_END":
                    log_reports.append(log_basename + "_aligned_PE_report.txt")
                    output_report = os.path.join("merged", sample.name, sample.name + ".merged_aligned_PE_report.txt")
                else:
                    log_reports.append(log_basename + "_aligned_SE_report.txt")
                    output_report = os.path.join("merged", sample.name, sample.name + ".merged_aligned_SE_report.txt")

            # Job creation
            mkdir_job = Job(command="mkdir -p merged/" + sample.name)
            merge_job = Job(log_reports, [output_report],
                            [['merge_bismark_alignment_report', 'module_python']],
                            command="""python {script_loc} -o {output} -n {name} {logs}""".format(
                                # A custom made script that should always be in the episeq directory, but not sure
                                # how to reference that directory. TODO: Find script automatically?
                                script_loc=os.path.abspath(config.param('DEFAULT', 'extern_script')),
                                output=output_report,
                                name=sample.name,
                                logs=' '.join(log_reports)))

            job = concat_jobs([mkdir_job, merge_job], name="merge_align_reports." + sample.name)
            jobs.append(job)
        return jobs

    def picard_merge_sam_files(self):
        """
        This step merges all readsets of each sample into one handy bam file. Here, if a readset is defined by a
        bam file, it will finally be used to add with other readsets. Because merging multiple alignments together
        can dramatically change the coverage content, it is recalculated to reflect the combined reads.

        While already defined in the bfx module, I want to use Picard Tools v2.0.1. This (and newer versions) has
        a different syntax than before, requiring me to rewrite the job definition.

        :return: A list of jobs that needs to be executed in this step.
        :rtype: list(Job)
        """

        jobs = []
        for sample in self.samples:
            # Bam files from pipeline (FASTQ)
            processed_fastq_pe = [os.path.join('aligned', sample.name, readset.name + "_aligned_pe.bam") for
                                  readset in sample.readsets if readset.run_type == 'PAIRED_END' and not readset.bam]
            processed_fastq_se = [os.path.join('aligned', sample.name, readset.name + "_aligned.bam") for
                                  readset in sample.readsets if readset.run_type == 'SINGLE_END' and not readset.bam]
            processed_fastq = processed_fastq_pe + processed_fastq_se
            # Bam files from user, if specified by readset file. Not exclusive with having fastq
            listed_bam_files = [readset.bam for readset in sample.readsets if readset.bam != '' and not readset.fastq1]

            input_files = processed_fastq + listed_bam_files  # All bam files that belong to the sample
            merge_prefix = os.path.join('merged', sample.name)
            output_bam = os.path.join(merge_prefix, sample.name + '.merged.bam')
            mkdir_job = Job(command='mkdir -p ' + merge_prefix)

            # I want to use Picard Tools v2.0.1, which has a different syntax than v1.x
            if len(input_files) > 1:
                picard_v2 = Job(
                    input_files,
                    [output_bam],
                    [
                        ['picard_merge_sam_files', 'module_java'],
                        ['picard_merge_sam_files', 'module_picard']
                    ],
                    command="""\
    java -Djava.io.tmpdir={tmp_dir} {java_other_options} -Xmx{ram} \\
      -jar $PICARD_HOME/picard.jar MergeSamFiles \\
      VALIDATION_STRINGENCY=SILENT \\
      TMP_DIR={tmp_dir} \\
      {inputs} \\
      OUTPUT={output} \\
      USE_THREADING=true \\
      SORT_ORDER=queryname \\
      MAX_RECORDS_IN_RAM={max_records_in_ram}""".format(
                        tmp_dir=config.param('picard_merge_sam_files', 'tmp_dir'),
                        java_other_options=config.param('picard_merge_sam_files', 'java_other_options'),
                        ram=config.param('picard_merge_sam_files', 'ram'),
                        inputs=" \\\n  ".join(["INPUT=" + in_put for in_put in input_files]),
                        output=output_bam,
                        max_records_in_ram=config.param('picard_merge_sam_files', 'max_records_in_ram', type='int')),
                    # removable_files=[output_bam, re.sub("\.([sb])am$", ".\\1ai", output_bam)],
                    local=config.param('picard_merge_sam_files', 'use_localhd', required=False))
                job = concat_jobs([mkdir_job, picard_v2], name="picard_merge_sam_files." + sample.name)

            elif len(input_files) == 1:  # Save time and resources by just copying the single data source
                input_nuc_stats = os.path.join('aligned', sample.name, sample.readsets[0].name)
                if processed_fastq_pe:
                    input_nuc_stats += '_aligned_PE_report.txt'
                else:
                    input_nuc_stats += '_aligned_SE_report.txt'
                target_readset_bam = input_files[0]
                job = concat_jobs([
                    mkdir_job,
                    Job([input_files[0]], [output_bam],
                        command="cp -s -L -f " + os.path.join(os.getcwd(), 'output', target_readset_bam) +
                                " " + os.path.join(os.getcwd(), 'output', output_bam),
                        removable_files=[output_bam])],
                    name="symlink_readset_sample_bam." + sample.name)
            else:
                raise ValueError('Sample ' + sample.name + ' has no readsets!')
            jobs.append(job)
        return jobs

    def merged_nuc_stats(self):
        """
        This step calculates the new nucleotide coverge for the merged bam file.
        This independent step is to reduce the likelihood of failing mid step.

        The output of this step is found in the relevant merge sub-directory.

        :return: A list of jobs that needs to be executed in this step.
        :rtype: list(Job)
        """
        jobs = []
        for sample in self.samples:
            output_dir = os.path.join('merged', sample.name)
            if os.path.exists(os.path.join(output_dir, sample.name) + '.merged.nucleotide_stats.txt'):
                continue
            merged_bam = os.path.join(output_dir, sample.name) + '.merged.bam'
            job = self.__bam2nuc__(output_dir, sample.name, '.merged', merged_bam)
            jobs.append(job)
        return jobs

    def bismark_deduplicate(self):
        """
        Calls the de-duplication module from Bismark. The module operates on the output alignment files from
        Bismark, creating a new output file (SAM default, BAM possible) with the suffix *.deduplicated.bam.
        The output file appears in the same directory as the input file, but this method will move the output file
        to its own directory.

        :return: A list of jobs that needs to be executed in this step.
        :rtype: list(Job)
        """

        jobs = []
        for sample in self.samples:
            work_dir = os.path.join('dedup', sample.name)
            in_file = os.path.join('merged', sample.name, sample.name + '.merged.bam')
            out_file = os.path.join(work_dir, sample.name + '.merged.deduplicated.bam')
            report_file = os.path.join(work_dir, sample.name + '.merged.deduplication_report.txt')
            protocol = sample.readsets[0].library
            run_type = sample.readsets[0].run_type

            # Get I/O
            if run_type == 'PAIRED_END':
                in_report_file = os.path.join('merged', sample.name, sample.name + '.merged_aligned_PE_report.txt')
                copy_report = os.path.join(work_dir, sample.name + '.merged.deduplication_aligned_PE_report.txt')
            else:
                in_report_file = os.path.join('merged', sample.name, sample.name + '.merged_aligned_SE_report.txt')
                copy_report = os.path.join(work_dir, sample.name + '.merged.deduplication_aligned_SE_report.txt')

            mkdir_job = Job(command='mkdir -p ' + work_dir)
            copy_job = Job([in_report_file],
                           [copy_report],
                           command="cp -fu " + in_report_file + " " + copy_report,
                           report_files=[copy_report])

            # A job name that is different from the heading will not use the params listed. (Uses default, instead)
            if protocol == 'RRBS':  # Deduplication is not recommended for RRBS datatypes. Keep what we have
                # You can only make a relative link in the current directory, so use absolute paths.
                abs_in_file = os.path.join(self.output_dir, in_file)
                abs_out_file = os.path.join(self.output_dir, out_file)
                job = concat_jobs([mkdir_job,
                                   Job([in_file], [out_file], command="cp -Ls -f " + abs_in_file + " " + abs_out_file),
                                   copy_job],
                                  name="skip_rrbs_deduplicate." + sample.name)
            else:
                merge_job = Job([in_file],
                                module_entries=[['bismark_deduplicate', 'module_samtools'],
                                                ['bismark_deduplicate', 'module_perl'],
                                                ['bismark_deduplicate', 'module_bismark']],
                                command="""deduplicate_bismark {type} --bam {other} {input}""".format(
                                    type='--paired' if run_type == 'PAIRED_END' else '--single',
                                    input=in_file,
                                    other=config.param('bismark_deduplicate', 'other_options', required=False)))
                move_bam = Job(output_files=[out_file],
                               command='mv -fu ' + os.path.join(os.path.dirname(in_file),
                                                                os.path.basename(out_file)) + ' ' + out_file,
                               removable_files=[out_file])
                move_log = Job(output_files=[report_file],
                               command='mv -fu ' + os.path.join(os.path.dirname(in_file),
                                                                os.path.basename(report_file)) + ' ' + report_file,
                               report_files=[report_file])
                job = concat_jobs([mkdir_job, merge_job, move_bam, move_log, copy_job],
                                  name='bismark_deduplicate.' + sample.name)
            jobs.append(job)
        return jobs

    def calc_dedup_nucleotide_coverage(self):
        """
        This step recalculates the nucleotide coverage values using Bismark's bam2nuc script. This recalculation is
        done after deduplication. The removal of some reads will likely alter this value, so this is helps update it.

        The output file is found in the related folder under dedup/.

        :return: A list of jobs that needs to be executed in this step.
        :rtype: list(Job)
        """
        jobs = []
        for sample in self.samples:
            if sample.readsets[0].library == 'RRBS':
                continue  # No deduplication done for RRBS samples
            output_dir = os.path.join('dedup', sample.name)
            merged_bam = os.path.join(output_dir, sample.name) + '.merged.deduplicated.bam'
            job = self.__bam2nuc__(output_dir, sample.name, '.merged.deduplicated', merged_bam)
            jobs.append(job)
        return jobs

    def bismark_methylation_caller(self):
        """
        This step extracts the methylation call for every single cytosine analyzed from the Bismark result files.
        The following input files are accepted:
            1.	Bismark result files from previous alignment step
            2.	BAM files (unsorted) from readset file

        Input: Merged sample files (merged/)
        Output: Methylation calls in BedGraph format. (methyl_calls/)

        :return: A list of jobs that needs to be executed in this step.
        :rtype: list(Job)
        """

        jobs = []

        for sample in self.samples:
            # Either select aligned sample from previous alignment step or aligned BAM/SAM files in readset file
            merged_sample = self.select_input_files([[readset.bam for readset in sample.readsets], [
                os.path.join("dedup", sample.name, sample.name + '.merged.deduplicated.bam')]])
            report_files = [
                os.path.join("methyl_calls", sample.name, sample.name + ".merged.deduplicated.bedGraph.gz"),
                os.path.join("methyl_calls", sample.name, sample.name + ".merged.deduplicated.bismark.cov.gz"),
                os.path.join("methyl_calls", sample.name, sample.name + ".merged.deduplicated.M-bias.txt"),
                os.path.join("methyl_calls", sample.name, sample.name + ".merged.deduplicated_splitting_report.txt"),
                os.path.join("methyl_calls", sample.name, sample.name + ".merged.deduplicated.CpG_report.txt.gz")]
            other_files = [
                os.path.join("methyl_calls", sample.name, "CHG_OB_" + sample.name + ".merged.deduplicated.txt.gz"),
                os.path.join("methyl_calls", sample.name, "CHG_OT_" + sample.name + ".merged.deduplicated.txt.gz"),
                os.path.join("methyl_calls", sample.name, "CHH_OB_" + sample.name + ".merged.deduplicated.txt.gz"),
                os.path.join("methyl_calls", sample.name, "CHH_OT_" + sample.name + ".merged.deduplicated.txt.gz"),
                os.path.join("methyl_calls", sample.name, "CpG_OT_" + sample.name + ".merged.deduplicated.txt.gz"),
                os.path.join("methyl_calls", sample.name, "CpG_OB_" + sample.name + ".merged.deduplicated.txt.gz")]
            run_type = sample.readsets[0].run_type
            job = Job(
                merged_sample + ['bismark_prepare_genome/'],
                report_files + other_files,
                [['bismark_methylation_caller', 'module_samtools'],
                 ['bismark_methylation_caller', 'module_perl'],
                 ['bismark_methylation_caller', 'module_bismark']],
                command="""\
mkdir -p {directory}
bismark_methylation_extractor {library_type} {other} --multicore {core} --output {directory} \
--bedGraph --cytosine_report --gzip --genome_folder {genome} {sample}
        """.format(
                    directory=os.path.join("methyl_calls", sample.name),
                    library_type="--paired-end" if run_type == "PAIRED_END" else "--single-end",
                    other=config.param("bismark_methylation_caller", "other_options"),
                    core=config.param('bismark_methylation_caller', 'cores'),
                    sample=" ".join(merged_sample),
                    genome=os.path.join(self.output_dir, 'bismark_prepare_genome')),
                report_files=report_files,
                removable_files=other_files,
                name="bismark_methylation_caller." + sample.name)
            jobs.append(job)
        return jobs

    def bismark_html_report_generator(self):
        """
        Generates the Bismark Report page by combining data from alignment, deduplication, methylation, and
        nucleotide coverage reports. The alignment report is a requirement while all others are optional.

        :return: A list of jobs that needs to be executed in this step.
        :rtype: list(Job)
        """
        # Report Generation Setup
        report_file = 'report/EpiSeq.bismark_html_report_generator.md'
        report_data = 'report/data/bismark_html_report_generator'
        template_string_file = 'bismark_summary_report/' + \
                               datetime.datetime.now().strftime('%Y_%m_%d') + "_template_var_strings.txt"
        fill_in_entry = '| {sample} | {report} | {methyl_calls} | {submission} | {completion} |'

        jobs = []
        module_list = [['bismark_html_report_generator', 'module_samtools'],
                       ['bismark_html_report_generator', 'module_perl'],
                       ['bismark_html_report_generator', 'module_bismark']]
        for sample in self.samples:
            report_list = ['', '', '', '', '']
            # Required file
            if sample.readsets[0].run_type == 'PAIRED_END':
                report_list[0] = os.path.join("merged", sample.name,
                                              sample.name + ".merged_aligned_PE_report.txt")
            else:
                report_list[0] = os.path.join("merged", sample.name,
                                              sample.name + ".merged_aligned_SE_report.txt")
            # Optional output depending on run type
            if sample.readsets[0].library != 'RRBS':
                report_list[1] = os.path.join('dedup', sample.name, sample.name + '.merged.deduplication_report.txt')
                report_list[4] = os.path.join('dedup', sample.name, sample.name +
                                              '.merged.deduplicated.nucleotide_stats.txt')
            else:
                report_list[1] = None
                report_list[4] = os.path.join('merged', sample.name, sample.name +
                                              '.merged.nucleotide_stats.txt')
            # Based on our pipeline, these files are always generated, so make it manditory
            report_list[2] = os.path.join("methyl_calls", sample.name, sample.name +
                                          ".merged.deduplicated_splitting_report.txt")
            report_list[3] = os.path.join("methyl_calls", sample.name, sample.name + ".merged.deduplicated.M-bias.txt")

            # Output file!
            html_report = os.path.join('bismark_summary_report', sample.name + '_final_bismark_report.html')

            # Job creation
            mkdir_job = Job(command='mkdir -p ' + os.path.dirname(html_report))
            job = Job(input_files=report_list, output_files=[html_report],
                      module_entries=module_list, report_files=[html_report],
                      command="""bismark2report -o {out} --verbose --alignment_report {align} \
                      {dedup} {split} {mbias} {nt}""".format(
                          out=html_report,
                          align=report_list[0],
                          dedup=' --dedup_report ' + report_list[1] if report_list[1] else '',
                          split=' --splitting_report ' + report_list[2] if report_list[2] else '',
                          mbias=' --mbias_report ' + report_list[3] if report_list[3] else '',
                          nt=' --nucleotide_report ' + report_list[4] if report_list[4] else ''))

            # Report Creation
            zip_file = os.path.join(report_data, sample.name + '_methyl.zip')
            report_entry = fill_in_entry.format(
                sample=sample.name,
                report='[View HTML](' + os.path.join(report_data, os.path.basename(html_report)) + ')',
                methyl_calls='[Download Raw Data](' + zip_file + ')',
                submission=datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
                completion="$(date +%Y-%m-%d %H:%M:%S)")
            command = """\
            # Mutual exclusion with other jobs
            flock -x ".{table_hold}.lock" -c "echo -e {entry} >> {table_hold}" && \\
            mkdir -p {data_dir} && \\
            cp -f {report} {data_dir} && \\
            zip {zip_file} {methyl_calls} && \\
            table=$(cat {table_hold}) && \\
            pandoc \\
            {report_template_dir}/{basename_report_file} \\
                --template {report_template_dir}/{basename_report_file} \\
                --variable data_table="$table" \\
                --to markdown \\
            > {report_file}
            """.format(
                table_hold=template_string_file,
                entry=report_entry,
                data_dir=report_data,
                methyl_calls='methyl_calls/' + sample.name + '/' + sample.name + '.merged.deduplicated.*',
                zip_file=zip_file,
                report=html_report,
                report_template_dir=self.report_template_dir,
                basename_report_file=os.path.basename(report_file),
                report_file=report_file
            )
            report_job = Job(output_files=[report_file,
                                           os.path.join(report_data, os.path.basename(html_report)),
                                           zip_file],
                             report_files=[report_file],
                             module_entries=[['bismark_html_report_generator', 'module_pandoc']],
                             command=command)
            jobs.append(concat_jobs([mkdir_job, job, report_job], name='bismark_report.' + sample.name))
        return jobs

    def differential_methylated_pos(self):
        """
        This step finds a list of differentially methylated CpG sites with respect to a categorical
        phenotype (controls vs. cases). The BedGraph files from the previous methylation calling step are first combined
        to a BSRaw object with the R package BiSeq. Then, the dmpFinder function from the R package minfi is used to
        compute a F-test statistic on the beta values for the assayed CpGs in each sample. A p-value is then returned
        for each site with the option of correcting them for multiple testing. Differential analysis is done for each
        contrast specified in the design file

        The values from the design files dictates how the samples are treated and compared.

        Input: Methylation data (methyl_calls/)
        Output: A CSV file in differential_methylated_positions/

        :return: A list of jobs that needs to be executed in this step.
        :rtype: list(Job)
        """

        jobs = []

        for contrast in self.contrasts:
            # Determine the control and case samples to include in the analysis from the contrast
            contrast_samples = [sample for sample in contrast.controls + contrast.treatments]
            cov_files = [os.path.join("methyl_calls", sample.name, sample.name + ".merged.deduplicated.bismark.cov.gz")
                         for sample in contrast_samples]
            sample_group = ["control" if sample in contrast.controls else "case" for sample in contrast_samples]
            dmps_file = os.path.join("differential_methylated_positions",
                                     contrast.name + "_RRBS_differential_methylated_pos.csv")
            if len(contrast.controls) == 0 or contrast.treatments == 0 or len(contrast_samples) <= 2:  # No 1v1 or less
                log.warn("Insufficient sample size to compare case and control. Skipping contrast: " + contrast.name)
                continue
            job = Job(
                cov_files,
                [dmps_file],
                [
                    ["differential_methylated_pos", "module_R"],
                    ["differential_methylated_pos", "module_mugqic_R_packages"]
                ],
                command="""\
mkdir -p {directory}
R --no-save --no-restore <<-EOF
suppressPackageStartupMessages(library(BiSeq))
suppressPackageStartupMessages(library(minfi))
rrbs <- readBismark(c{samples}, colData=DataFrame(group=factor(c{group}), row.names=c{sample_names}))
rrbs.filtered <- rrbs[apply(totalReads(rrbs), 1, function(x) any(x > {coverage})),]
beta <- methLevel(rawToRel(rrbs.filtered))

#Use M values to do statistical tests because they are more reliable
#dmpFinder does not work with M values that are 0 or INF so the beta values must be shifted slightly
#Although there is no such thing as a beta value > 1, it will not matter in this step because only
#the average beta values are shown to the user
beta[beta == 0] = 0.000001
beta[beta == 1] = 0.999999
M <- log2(beta/(1-beta))

dmp <- dmpFinder(M, pheno=colData(rrbs.filtered)[,"group"], type="categorical")
dmp["pval"] <- p.adjust(dmp[,"pval"], method = "{padjust_method}")
dmp <- dmp[dmp["pval"] < {pvalue},]["pval"]

controls <- c({controls})
cases <- c({cases})
result = as.data.frame(rowRanges(rrbs.filtered))[1:4]
result["Avg Control Beta"] = rowMeans(beta[,controls])
result["Avg Case Beta"] = rowMeans(beta[,cases])
result["Avg Delta Beta"] = result[,"Avg Case Beta"] - result[,"Avg Control Beta"]
result <- merge(result, dmp, by=0)
result <- result[abs(result["Avg Delta Beta"]) > {delta_beta_threshold}]

write.csv(result, file="{dmps_file}", quote=FALSE, row.names=FALSE)

EOF
zip {zip_file} {dmps_file}
                """.format(
                    directory=os.path.dirname(dmps_file),
                    samples=tuple(cov_files),
                    group=tuple(sample_group),
                    sample_names=tuple([sample.name for sample in contrast_samples]),
                    coverage=config.param("differential_methylated_pos", "read_coverage"),
                    controls=', '.join(["'" + sample.name + "'" for sample in contrast.controls]),
                    cases=', '.join(["'" + sample.name + "'" for sample in contrast.treatments]),
                    padjust_method=config.param("differential_methylated_pos", "padjust_method"),
                    pvalue=config.param("differential_methylated_pos", "pvalue", type="float"),
                    delta_beta_threshold=config.param("differential_methylated_pos", "delta_beta_threshold",
                                                      type="float"),
                    dmps_file=dmps_file,
                    zip_file='report/data/differential_methylated_pos.zip'
                ),
                name="differential_methylated_pos." + contrast.name)

            jobs.append(job)

        return jobs

    def differential_methylated_regions(self):
        """
        Similar to differential_methylated_positions, this step looks at methylation patterns on a larger, regional
        level. This step compares large-scale differences in methylation as opposed to comparing local methylation
        sites.

        Input: Methylation data (methyl_calls/)
        Output: A CSV file in differential_methylated_regions/

        :return: A list of jobs that needs to be executed in this step.
        :rtype: list(Job)
        """
        jobs = []

        for contrast in self.contrasts:
            # Determine the control and case samples to include in the analysis from the contrast
            contrast_samples = [sample for sample in contrast.controls + contrast.treatments]
            cov_files = [os.path.join("methyl_calls", sample.name, sample.name +
                                      ".merged.deduplicated.bismark.cov.gz") for sample in contrast_samples]
            sample_group = ["control" if sample in contrast.controls else "case" for sample in contrast_samples]
            dmrs_file = os.path.join("differential_methylated_regions",
                                     contrast.name + "_RRBS_differential_methylated_regions.csv")

            job = Job(
                cov_files,
                [dmrs_file],
                [
                    ["differential_methylated_regions", "module_R"],
                    ["differential_methylated_regions", "module_mugqic_R_packages"]
                ],
                command="""\
mkdir -p {directory}
R --no-save --no-restore <<-EOF
suppressPackageStartupMessages(library(bumphunter))
suppressPackageStartupMessages(library(BiSeq))
library(doParallel)
registerDoParallel(cores={cores})

rrbs <- readBismark(c{samples}, colData=DataFrame(group=c{group}, row.names=c{sample_names}))
rrbs.filtered <- rrbs[apply(totalReads(rrbs), 1, function(x) any(x > {coverage})),]
beta <- methLevel(rawToRel(rrbs.filtered))
chr <- as.character(seqnames(rowRanges(rrbs.filtered)))
pos <- start(ranges(rowRanges(rrbs.filtered)))
pheno <- colData(rrbs.filtered)[,"group"]
designM <- model.matrix(~pheno)

dmrs <- bumphunterEngine(beta,
                         chr=chr,
                         pos=pos,
                         design=designM,
                         cutoff={delta_beta_threshold},
                         pickCutoffQ=0.99,
                         null_method=c("permutation","bootstrap"),
                         smooth=FALSE,
                         smoothFunction=locfitByCluster,
                         B={permutations},
                         verbose=TRUE,
                         maxGap=500)

dmrs <- na.omit(dmrs)

write.csv(dmrs\$table, "{dmrs_file}", quote=FALSE, row.names=FALSE)

EOF

zip {zip_file} {dmrs_file}
                """.format(
                    directory=os.path.dirname(dmrs_file),
                    samples=tuple(cov_files),
                    group=tuple(sample_group),
                    coverage=config.param("differential_methylated_regions", "read_coverage", type="int"),
                    sample_names=tuple([sample.name for sample in contrast_samples]),
                    cores=config.param('bismark_methylation_caller', 'cluster_cpu').split('=')[-1],
                    delta_beta_threshold=config.param("differential_methylated_regions", "delta_beta_threshold",
                                                      type="float"),
                    permutations=config.param("differential_methylated_regions", "permutations", type="int"),
                    dmrs_file=dmrs_file,
                    zip_file='report/data/differential_methylated_regions.zip'
                ),
                name="differential_methylated_regions." + contrast.name)

            jobs.append(job)

        return jobs

    @staticmethod
    def __bam2nuc__(output_dir, sample_name, suffix, in_bam):
        """
        Generates jobs for Bismark's bam2nuc script.

        :param output_dir: A specified output directory for the report file
        :type output_dir: str
        :param sample_name: The sample of the sample to run
        :type sample_name: str
        :param suffix: A suffix to add to the filename before '.nucleotide_stats.txt'
        :type suffix: str
        :param in_bam: The bam file to analyse.
        :type in_bam: str
        :return: A Job object to run bam2nuc
        :rtype: Job
        """
        output_file = os.path.join(output_dir, sample_name + suffix + '.nucleotide_stats.txt')
        coverage_calc = Job(
            input_files=[in_bam],
            output_files=[output_file],
            module_entries=[['bismark_deduplicate', 'module_samtools'],
                            ['bismark_deduplicate', 'module_perl'],
                            ['bismark_deduplicate', 'module_bismark']],
            command='bam2nuc --dir ' + output_dir + ' --genome_folder bismark_prepare_genome ' + in_bam,
            name='bam2nuc.' + sample_name)
        return coverage_calc

    @property
    def steps(self):
        """
        These are the steps to the pipeline and are listed in the order that should be processed

        :return: list
        """
        return [
            self.bismark_prepare_genome,
            self.pre_qc_check,
            self.trim_galore,
            self.bismark_align,
            self.merge_bismark_alignment_report,
            self.picard_merge_sam_files,
            self.merged_nuc_stats,
            self.bismark_deduplicate,
            self.calc_dedup_nucleotide_coverage,
            self.bismark_methylation_caller,
            self.bismark_html_report_generator,
            self.differential_methylated_pos,
            self.differential_methylated_regions
        ]


if __name__ == '__main__':
    Episeq()
