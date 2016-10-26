#!/usr/bin/env python

# Epigenetics pipeline for RRBS/WGSB data

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
    The Episeq pipeline takes FASTQ or BAM files (unsorted) as input
    """

    def __init__(self):
        self.argparser.add_argument("-d", "--design", help="design file", type=file)
        super(Episeq, self).__init__()

    def merge_fastq(self):

        """
        This step merges FASTQ files for samples with multiple readsets or creates symbolic links of FASTQ
        files for samples with only a single readset
        """

        jobs = []
        for sample in self.samples:
            merge_directory = os.path.join("merged", sample.name)
            merge_prefix = os.path.join(merge_directory, sample.name)
            mkdir_job = Job(command="mkdir -p " + merge_directory)
            # If one readset in a sample is of a particular type,
            # then all readsets are assumed to be of that type as well
            run_type = sample.readsets[0].run_type

            job = []
            # Samples with only one readset do not require to be merged.
            # A symbolic link to the FASTQ files are made instead
            if len(sample.readsets) == 1:
                target_fastq1 = sample.readsets[0].fastq1
                if run_type == "PAIRED_END":
                    target_fastq2 = sample.readsets[0].fastq2
                elif run_type == "SINGLE_END":
                    target_fastq2 = ""

                job = concat_jobs([
                    mkdir_job,
                    Job([target_fastq1], [merge_prefix + "_R1.fastq.gz"],
                        command="ln -s -f " + target_fastq1 + " " + merge_prefix + "_R1.fastq.gz"),
                    Job([target_fastq2], [merge_prefix + "_R2.fastq.gz" if run_type == "PAIRED_END" else ""],
                        command="ln -s -f " + target_fastq2 + " " + merge_prefix + "_R2.fastq.gz")
                ], name="symlink_fastq." + sample.name)

            # Samples with more than one readsets will be merged by their forward and/or reverse FASTQ files
            elif len(sample.readsets) > 1:
                fastq1 = [readset.fastq1 for readset in sample.readsets]
                if run_type == "PAIRED_END":
                    fastq2 = [readset.fastq2 for readset in sample.readsets]
                elif run_type == "SINGLE_END":
                    fastq2 = ""

                job = concat_jobs([
                    mkdir_job,
                    Job([fastq for fastq in fastq1],
                        [merge_prefix + "_R1.fastq.gz"],
                        command="cat {fastq1} > {merged_fastq1}".format(
                            fastq1=" ".join(fastq1),
                            merged_fastq1=merge_prefix + "_R1.fastq.gz")),
                    Job([fastq for fastq in fastq2],
                        [merge_prefix + "_R2.fastq.gz" if run_type == "PAIRED_END" else ""],
                        command="cat {fastq2} > {merged_fastq2}".format(
                            fastq2=" ".join(fastq2),
                            merged_fastq2=merge_prefix + "_R2.fastq.gz"))
                ], name="merge_fastq." + sample.name)

            jobs.append(job)

        return jobs

    def bismark_prepare_genome(self):
        """
        Bismark requires a processed reference genome to compare with the epigenome.
        """
        ref_seq = config.param('bismark_prepare_genome', 'genome_file', type='filepath')
        local_ref_seq = os.path.join('bismark_prepare_genome', os.path.basename(ref_seq))
        output_idx = "Bisulfite_Genome"

        if not os.path.isfile(local_ref_seq):
            run_job = Job([ref_seq], [output_idx],
                          module_entries=[["bismark_prepare_genome", "module_bowtie2"],
                                          ["bismark_prepare_genome", "module_samtools"]],
                          command="""\
    cp {src} .
    module load bismark/0.15
    bismark_genome_preparation {verbose} --bowtie2 .""".format(
                              src=ref_seq,
                              verbose='--verbose' if config.param('bismark_genome_preparation', 'verbose_logging',
                                                                  required=False, type='boolean') else ''),
                          name="bismark_prepare_genome")
        else:
            run_job = Job([ref_seq], [output_idx],
                          module_entries=[["bismark_prepare_genome", "module_bowtie2"],
                                          ["bismark_prepare_genome", "module_samtools"]],
                          command="""\
    module load bismark/0.15
    bismark_genome_preparation {verbose} --bowtie2 .""".format(
                              verbose='--verbose' if config.param('bismark_genome_preparation', 'verbose_logging',
                                                                  required=False, type='boolean') else ''),
                          name="bismark_prepare_genome")

        return [run_job]

    def trim_galore(self):

        """
        This step trims raw FASTQ files for quality control using Trim Galore! and runs fastqc at the end,
        courtesy of trim galore.
        """

        jobs = []
        for sample in self.samples:
            for readset in sample.readsets:
                run_type = readset.run_type
                protocol = readset.library
                trim_directory = os.path.join("trimmed", sample.name, readset.name)
                fq1_out = os.path.join(trim_directory, readset.name)
                output_files = []
                report_logs = []

                # Trim Galoree has no built in option to change the filenames of the output
                # Below are the default output names when running in paired or single mode
                if run_type == "PAIRED_END":
                    input_files = [readset.fastq1, readset.fastq2]
                    output_files = [fq1_out + "_1_val_1.fq.gz", fq1_out + "_2_val_2.fq.gz"]
                    report_logs = [fq1_out + '_1_trimming_report.txt', fq1_out + "_1_val_1_fastqc.html",
                                   fq1_out + '_2_trimming_report.txt', fq1_out + "_2_val_2_fastqc.html"]
                elif run_type == "SINGLE_END":
                    input_files = [readset.fastq1]
                    output_files = [fq1_out + "_trimmed.fq.gz"]
                    report_logs = [fq1_out + '_trimming_report.txt', fq1_out + "_fastqc.html"]

                mkdir_job = Job(command="mkdir -p " + trim_directory)
                job = concat_jobs([
                    mkdir_job,
                    Job(
                        input_files,
                        output_files,
                        [['trim_galore', 'module_fastqc']],
                        command="""\
    module load trim_galore/0.4.1
    module load cutadapt/1.10
    trim_galore {protocol} {library_type} {non_directional} {other} --output_dir {directory} \
    --fastqc_args "-t 4" {fastq1} {fastq2}
        """.format(
                            library_type="--paired" if run_type == "PAIRED_END" else "",
                            protocol='--rrbs' if protocol == 'RRBS' else '',
                            non_directional='--non_directional' if run_type == 'PAIRED_END' and protocol == 'RRBS'
                            else '',
                            other=config.param("trim_galore", "other_options"),
                            directory=trim_directory,
                            fastq1=input_files[0],
                            fastq2='' if run_type == "SINGLE_END" else input_files[1]
                        ),
                        report_files=report_logs
                    )], name="trim_galore." + readset.name)
                jobs.append(job)
        return jobs

    def bismark_align(self):
        """
        This step aligns trimmed reads to a bisulfite converted reference genome using Bismark [link]
        """

        # Multicore funtionality for Bismark is currently not supported when using the option --basename
        jobs = []
        for sample in self.samples:
            for readset in sample.readsets:
                trim_prefix = os.path.join("trimmed", sample.name, readset.name)
                align_directory = os.path.join("aligned", sample.name)
                readset_sam = os.path.join(align_directory, readset.name + "_aligned_pe.bam")
                report_log = [os.path.join(align_directory, readset.name + "aligned_PE_report.txt")]
                run_type = readset.run_type
                protocol = readset.library

                if run_type == "PAIRED_END":
                    input_files = [os.path.join(trim_prefix, readset.name + "_1_val_1.fq.gz"),
                                   os.path.join(trim_prefix, readset.name + "_2_val_2.fq.gz")]
                elif run_type == "SINGLE_END":
                    input_files = [os.path.join(trim_prefix, readset.name + "_trimmed.fq.gz")]

                mkdir_job = Job(command="mkdir -p " + align_directory)
                job = concat_jobs([
                    mkdir_job,
                    Job(
                        input_files + ["Bisulfite_Genome"],
                        [readset_sam],
                        [["bismark_align", "module_bowtie2"],
                         ["bismark_align", "module_samtools"]],
                        # ['bismark_align', 'module_perl']],
                        # Do not import. Bismark is expecting /usr/bin/perl. Will raise errors otherwise.
                        command="""\
module load bismark/0.15
bismark -q {directional} {other} --output_dir {directory} --basename {basename} --genome_folder . {fastq1} {fastq2}
        """.format(
                            directory=align_directory,
                            other=config.param("bismark_align", "other_options"),
                            fastq1=input_files[0],
                            fastq2=input_files[1] if run_type == "PAIRED_END" else "--single_end",
                            directional='--non_directional' if protocol == 'RRBS' else '',
                            basename=readset.name + '_aligned'
                        ),
                        report_files=report_log
                    )], name="bismark_align." + readset.name)
                jobs.append(job)
        return jobs

    def merge_bismark_alignment_report(self):
        """

        :return:
        :rtype:
        """
        search_total_seqs = '^Sequence(?:s|pairs) analysed in total:\s+([0-9]+)'
        search_uniq_hit = '^Number of (?:paired-end)? alignments with a unique best hit:\s+([0-9]+)'
        search_no_align = 'with no alignments under any condition:\s+([0-9]+)'
        search_not_uniq = 'did not map uniquely:\s+([0-9]+)'
        search_discarded = 'were discarded because genomic sequence.*?:\s+([0-9]+)'
        search_ct_ct = 'CT/GA/CT:\s+([0-9]+)\s+.*'
        search_ga_ct = 'GA/CT/CT:\s+([0-9]+)\s+.*'
        search_ga_ga = 'GA/CT/GA:\s+([0-9]+)\s+.*'
        search_ct_ga = 'CT/GA/GA:\s+([0-9]+)\s+.*'
        search_direct = 'complementary strands being rejected in total:\s+([0-9]+)'
        search_totalC = 'Total number of C.*?\s+([0-9]+)'
        search_mC_cpg = 'methylated (?:C\'s)? in CpG context:\s+([0-9]+)'
        search_mC_chg = 'methylated (?:C\'s)? in CHG context:\s+([0-9]+)'
        search_mC_chh = 'methylated (?:C\'s)? in CHH context:\s+([0-9]+)'
        search_mC_unk = 'methylated (?:C\'s)? in Unknown context:\s+([0-9]+)'
        search_C_cpg = '(?:unmethylated C\'s)|(?:C to T conversions) in CpG context:\s+([0-9]+)'
        search_C_chg = '(?:unmethylated C\'s)|(?:C to T conversions) in CHG context:\s+([0-9]+)'
        search_C_chh = '(?:unmethylated C\'s)|(?:C to T conversions) in CHH context:\s+([0-9]+)'
        search_C_unk = '(?:unmethylated C\'s)|(?:C to T conversions) in Unknown context:\s+([0-9]+)'
        search_all = [search_total_seqs, search_totalC, search_direct,
                      search_uniq_hit, search_no_align, search_not_uniq, search_discarded,
                      search_ct_ct, search_ga_ct, search_ga_ga, search_ct_ga,
                      search_mC_cpg, search_mC_chg, search_mC_chh, search_mC_unk,
                      search_C_cpg, search_C_chg, search_C_chh, search_C_unk]

        for sample in self.samples:
            log_reports = []
            output_dir = os.path.join("reports", sample.name)
            output_report = os.path.join(output_dir, sample.name + "_aligned_report.txt")
            align_directory = os.path.join("aligned", sample.name)

            # Get all log reports for this sample
            for readset in sample.readsets:
                log_basename = os.path.join(align_directory, readset.name)
                if readset.run_type == "PAIRED_END":
                    log_reports.append(log_basename + "_aligned_PE_report.txt")
                else:
                    log_reports.append(log_basename + "_aligned_SE_report.txt")

            # Set variables
            total_seqs, total_C, direction_rejected = (0, 0, 0)
            uniq_hit, no_align, not_uniq, discarded = (0, 0, 0, 0)
            ct_ct, ga_ct, ga_ga, ct_ga = (0, 0, 0, 0)
            mC_cpg, mC_chg, mC_chh, mC_unk = (0, 0, 0, 0)
            C_cpg, C_chg, C_chh, C_unk = (0, 0, 0, 0)
            value_all = [total_seqs, total_C, direction_rejected,
                         uniq_hit, no_align, not_uniq, discarded,
                         ct_ct, ga_ct, ga_ga, ct_ga,
                         mC_cpg, mC_chg, mC_chh, mC_unk,
                         C_cpg, C_chg, C_chh, C_unk]

            with [open(logs).readlines() for logs in log_reports] as file_data:
                for each_file in file_data:
                    for line in each_file:
                        for (val, total) in zip(search_all, value_all):
                            result = re.search(val, line)
                            if result:
                                total += result.group(1)
                                continue
            with open(output_report, 'w') as writer:
                writer.write("""
Merged Bismark report for: {readsets}

Final Alignment report
======================

Sequences analysed in total:\t{seqs}
Number of alignments with a unique hit:\t{uniq}
Mapping efficiency:\t{efficient}
Sequences with no alignments under any condition:\t{no_align}
Sequences did not map uniquely:\t{no_uniq}
Sequences which were discarded because genomic sequence could not be extracted:\t{discarded}

Number of sequences with unique best (first) alignment came from the bowtie output:
CT/GA/CT:\t{ct_ct}\t((converted) top strand)
GA/CT/CT:\t{ga_ct}\t(complementary to (converted) top strand)
GA/CT/GA:\t{ga_ga}\t(complementary to (converted) bottom strand)
CT/GA/GA:\t{ct_ga}\t((converted) bottom strand)

Number of alignments to (merely theoretical) complementary strands being rejected in total:\t{reject}

Final Cytosine Methylation Report
=================================
Total number of C's analysed:\t{total_c}

Total methylated C's in CpG context:\t{mc_cpg}
Total methylated C's in CHG context:\t{mc_chg}
Total methylated C's in CHH context:\t{mc_chh}
Total methylated C's in Unknown context:\t{mc_unk}


Total unmethylated C's in CpG context:\t{c_cpg}
Total unmethylated C's in CHG context:\t{c_chg}
Total unmethylated C's in CHH context:\t{c_chh}
Total unmethylated C's in Unknown context:\t{c_unk}

""".format(readsets=' '.join(sample.readsets), seqs=total_seqs,
           uniq=uniq_hit, efficient='', no_align=no_align, no_uniq=not_uniq, discarded=discarded,
           ct_ct=ct_ct, ga_ct=ga_ct, ga_ga=ga_ga, ct_ga=ct_ga, reject=direction_rejected,
           total_c=total_C, mc_cpg=mC_cpg, mc_chg=mC_chg, mc_chh=mC_chh, mc_unk=mC_unk,
           c_cpg=C_cpg, c_chg=C_chg, c_chh=C_chh, c_unk=C_unk))

    def picard_merge_sam_files(self):
        """

        :return:
        :rtype:
        """

        jobs = []
        for sample in self.samples:
            readsets = [os.path.join('aligned', sample.name,
                                     readset.name + "_aligned_pe.bam") for readset in sample.readsets]
            merge_prefix = 'merged'
            output_bam = os.path.join(merge_prefix, sample.name + '.merged.bam')

            mkdir_job = Job(command='mkdir -p ' + merge_prefix)

            # I want to use Picard Tools v2.0.1, which has a different syntax than v1.x
            if len(sample.readsets) > 1:
                picard_v2 = Job(
                    readsets,
                    [output_bam, re.sub("\.([sb])am$", ".\\1ai", output_bam)],
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
                        inputs=" \\\n  ".join(["INPUT=" + in_put for in_put in readsets]),
                        output=output_bam,
                        max_records_in_ram=config.param('picard_merge_sam_files', 'max_records_in_ram', type='int')),
                    removable_files=[output_bam, re.sub("\.([sb])am$", ".\\1ai", output_bam)],
                    local=config.param('picard_merge_sam_files', 'use_localhd', required=False))
                job = concat_jobs([mkdir_job, picard_v2],
                                  name="picard_merge_sam_files." + sample.name)  # Name must be set to match picard

            elif len(sample.readsets) == 1:
                readset_bam = readsets[0]
                if os.path.isabs(readset_bam):
                    target_readset_bam = readset_bam
                else:
                    target_readset_bam = os.path.relpath(readset_bam, merge_prefix)
                readset_index = re.sub("\.bam$", ".bai", readset_bam)
                target_readset_index = re.sub("\.bam$", ".bai", target_readset_bam)
                output_idx = re.sub("\.bam$", ".bai", output_bam)

                job = concat_jobs([
                    mkdir_job,
                    Job([readset_bam], [output_bam], command="ln -s -f " + target_readset_bam + " " + output_bam,
                        removable_files=[output_bam]),
                    Job([readset_index], [output_idx],
                        command="ln -s -f " + target_readset_index + " " + output_idx, removable_files=[output_idx])],
                    name="symlink_readset_sample_bam." + sample.name)
            else:
                raise ValueError('Sample ' + sample.name + ' has no readsets!')
            jobs.append(job)
        return jobs

    def bismark_methylation_caller(self):
        """
        This step extracts the methylation call for every single cytosine analyzed from the Bismark result files.
        The following input files are accepted:
            1.	Bismark result files from previous alignment step
            2.	BAM files (unsorted) from readset file
        """

        jobs = []

        for sample in self.samples:
            # Either select aligned sample from previous alignment step or aligned BAM/SAM files in readset file
            merged_sample = self.select_input_files([[readset.bam for readset in sample.readsets], [
                os.path.join("merged", sample.name + ".merged.bam")]])
            run_type = sample.readsets[0].run_type
            job = Job(
                merged_sample,
                [os.path.join("methyl_calls", sample.name, sample.name + "_aligned_pe.sam.bismark.cov.gz")],
                [['bismark_methylation_caller', 'module_samtools']],
                command="""\
        mkdir -p {directory}
        module load bismark/0.15
        bismark_methylation_extractor {library_type} {other} --output {directory} --bedGraph {sample}
        """.format(
                    directory="methyl_calls",
                    library_type="--paired-end" if run_type == "PAIRED_END" else "--single-end",
                    other=config.param("bismark_methylation_caller", "other_options"),
                    sample=" ".join(merged_sample)),
                name="bismark_methylation_caller." + sample.name)

            jobs.append(job)

        return jobs

    def differential_methylated_pos(self):

        """
        This step finds a list of differentially methylated CpG sites with respect to a categorical
        phenotype (controls vs. cases). The BedGraph files from the previous methylation calling step are first combined
        to a BSRaw object with the R package BiSeq. Then, the dmpFinder function from the R package minfi is used to
        compute a F-test statistic on the beta values for the assayed CpGs in each sample. A p-value is then returned
        for each site with the option of correcting them for multiple testing. Differential analysis is done for each
        contrast specified in the design file
        """

        jobs = []

        for contrast in self.contrasts:
            # Determine the control and case samples to include in the analysis from the contrast
            contrast_samples = [sample for sample in contrast.controls + contrast.treatments]
            cov_files = [os.path.join("methyl_calls", sample.name, sample.name + "_aligned_pe.sam.bismark.cov.gz") for
                         sample in contrast_samples]
            sample_group = ["control" if sample in contrast.controls else "case" for sample in contrast_samples]
            dmps_file = os.path.join("differential_methylated_positions",
                                     contrast.name + "_RRBS_differential_methylated_pos.csv")

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
        beta[beta == 0] = 0.0001
        beta[beta == 1] = 1.0001
        M <- log2(beta/(1-beta))

        dmp <- dmpFinder(M, pheno=colData(rrbs.filtered)[,"group"], type="categorical")
        dmp["pval"] <- p.adjust(dmp[,"pval"], method = "{padjust_method}")
        dmp <- dmp[dmp["pval"] < {pvalue},]["pval"]

        controls <- c{controls}
        cases <- c{cases}
        result = as.data.frame(rowRanges(rrbs.filtered))[1:4]
        result["Avg Control Beta"] = rowMeans(beta[,controls])
        result["Avg Case Beta"] = rowMeans(beta[,cases])
        result["Avg Delta Beta"] = result[,"Avg Case Beta"] - result[,"Avg Control Beta"]
        result <- merge(result, dmp, by=0)
        result <- result[abs(result["Avg Delta Beta"]) > {delta_beta_threshold}]

        write.csv(result, file="{dmps_file}", quote=FALSE, row.names=FALSE)

        EOF""".format(
                    directory=os.path.dirname(dmps_file),
                    samples=tuple(cov_files),
                    group=tuple(sample_group),
                    sample_names=tuple([sample.name for sample in contrast_samples]),
                    coverage=config.param("differential_methylated_pos", "read_coverage"),
                    controls=tuple([sample.name for sample in contrast.controls]),
                    cases=tuple([sample.name for sample in contrast.treatments]),
                    padjust_method=config.param("differential_methylated_pos", "padjust_method"),
                    pvalue=config.param("differential_methylated_pos", "pvalue", type="float"),
                    delta_beta_threshold=config.param("differential_methylated_pos", "delta_beta_threshold",
                                                      type="float"),
                    dmps_file=dmps_file
                ),
                name="differential_methylated_pos." + contrast.name)

            jobs.append(job)

        return jobs

    def differential_methylated_regions(self):

        """

        :return:
        :rtype:
        """
        jobs = []

        for contrast in self.contrasts:
            # Determine the control and case samples to include in the analysis from the contrast
            contrast_samples = [sample for sample in contrast.controls + contrast.treatments]
            cov_files = [os.path.join("methyl_calls", sample.name + "_aligned_pe.sam.bismark.cov.gz") for sample in
                         contrast_samples]
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

        EOF""".format(
                    directory=os.path.dirname(dmrs_file),
                    samples=tuple(cov_files),
                    group=tuple(sample_group),
                    coverage=config.param("differential_methylated_regions", "read_coverage", type="int"),
                    sample_names=tuple([sample.name for sample in contrast_samples]),
                    cores=config.param('bismark_methylation_caller', 'cluster_cpu').split('=')[-1],
                    delta_beta_threshold=config.param("differential_methylated_regions", "delta_beta_threshold",
                                                      type="float"),
                    permutations=config.param("differential_methylated_regions", "permutations", type="int"),
                    dmrs_file=dmrs_file
                ),
                name="differential_methylated_regions." + contrast.name)

            jobs.append(job)

        return jobs

    @property
    def steps(self):
        """

        :return: list
        """
        return [
            self.bismark_prepare_genome,
            self.trim_galore,
            self.bismark_align,
            self.picard_merge_sam_files,
            self.bismark_methylation_caller,
            self.differential_methylated_pos,
            self.differential_methylated_regions
        ]


if __name__ == '__main__':
    Episeq()
