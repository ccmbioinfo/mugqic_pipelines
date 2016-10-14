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


def add_pg_file(pgFile, cmdID, cmdVer, cmdLine, cmdPP):
    
    return Job(
        [None],
        [pgFile],
        [
            ["add_pg_file", "module_perl"],
            ["add_pg_file", "module_mugqic_tools"]
        ],
        command="""\
{script_path}/scripts/addToPGFile.pl \\
  --cmdID {cmdID} \\
  --cmdVer {cmdVer} \\
  --cmdLine {cmdLine} \\
  --cmdPP {cmdPP} \\
  --pgFile {pgFile}""".format(
            script_path=config.param("DEFAULT", "forge_location"),
            cmdID=cmdID,
            cmdVer=cmdVer,
            cmdLine=cmdLine,
            cmdPP=cmdPP,
            pgFile=pgFile
    ))


def add_to_vcf_header(input, output, header_command):
    
    samtools_ver = config.param("DEFAULT", "module_samtools").split("/")[-1]
    samtools_options = config.param("samtools_mpileup", "mpileup_other_options")
    samtools_ref = config.param("samtools_mpileup", "genome_fasta", type="filepath")

    new_header="##PG:samtools_mpileup,VN:" + samtools_ver + ",CL:" + header_command + " > " + output

    return Job(
        [None],
        [output],
        [
            ["add_to_vcf_header", "module_perl"],
            ["add_to_vcf_header", "module_mugqic_tools"]
        ],
        command="""\
{script_path}/scripts/add_to_vcf_header.pl \\
  --add "{new_header}" \\
  > {output}""".format(
            script_path=config.param("DEFAULT", "forge_location"),
            new_header=new_header,
            output=output
    ))
            

def aggregate_coverage(input, output, other_options=None):
    
    ccds_genes = config.param("gatk_depth_of_coverage", "ccds_genes_for_gatk")

    return Job(
        [input],
        [output],
        [
            ["aggregate_coverage", "module_perl"],
            ["aggregate_coverage", "module_mugqic_tools"]
        ],
        command="""\
{script_path}/scripts/aggregateCoverageByGene.pl {other_options}\\
  --ccdsGenes {ccds_genes} \\
  --sample {input} \\
  > {output}""".format(
            script_path=config.param("DEFAULT", "forge_location"),
            input=input,
            output=output,
            other_options=other_options + " \\\n  " if other_options else "",
            ccds_genes=ccds_genes
    ))


def allele_ratio_hist(input, pdf, output):

    return Job(
        [input],
        [output],
        [
            ["allele_ratio_hist", "module_perl"],
            ["allele_ratio_hist", "module_R"],
            ["allele_ratio_hist", "module_mugqic_tools"]
        ],
        command="""\
{script_path}/scripts/MakeHistAlleleRatio \\
  --input {input} \\
  --pdf {pdf} > {output}""".format(
            script_path=config.param("DEFAULT", "forge_location"),
            input=input,
            pdf=pdf,
            output=output
        ))


def allele_ratio_metrics(input, output, ini_section="allele_ratio_metrics"):

    return Job(
        [input],
        [output],
        [
            [ini_section, "module_perl"],
            [ini_section, "module_R"],
            [ini_section, "module_mugqic_tools"]
        ],
        command="""\
{script_path}/scripts/allele_ratio_metrics.pl \\
  --vcf {input} \\
  --minReadCount {minReadCount} \\
  --minQ {minQuality} \\
  --minMapQ {minMapQ} \\
  --prevSeenThreshold {prevSeenThresh} \\
  --maxMAF {maxMAF} \\
  --filterErrors > {output}""".format(
            script_path=config.param("DEFAULT", "forge_location"),
            input=input,
            minReadCount=config.param(ini_section, "MIN_READ_COUNT"),
            minQuality=config.param(ini_section, "MIN_QUALITY"),
            minMapQ=config.param(ini_section, "MIN_MAPQ"),
            prevSeenThresh=config.param(ini_section, "REMOVE_PREV_SEEN"),
            maxMAF=config.param(ini_section, "MAFThreshold"),
            output=output
        ))



def evs(input, output, ini_section="evs"):

    return Job(
        [input],
        [output],
        [
            [ini_section, "module_perl"],
            [ini_section, "module_mugqic_tools"]
        ],
        command="""\
{script_path}/scripts/annotate_evs.pl \\
  --vcf {input} \\
  --evsdb {evsdb} > {output} """.format(
            script_path=config.param("DEFAULT", "forge_location"),
            evsdb=config.param(ini_section, "evsdb", type="filepath"),
            input=input,
            output=output
    ))


def filter(input, output, options):
    
    return Job(
        [input],
        [output],
        [
            ["filter", "module_perl"],
            ["filter", "module_mugqic_tools"]
        ],
        command="""\
{script_path}/scripts/filter_combined_variants.pl \\
  --vcf {input}{options} > {output}""".format(
            script_path=config.param("DEFAULT", "forge_location"),
            input=input,
            options=options,
            output=output
        ))


def gene_mutation_counts(input, output, ini_section="gene_mutation_counts"):
    
    return Job(
        [input],
        [output],
        [
            [ini_section, "module_perl"],
            [ini_section, "module_mugqic_tools"]
        ],
        command="""\
{script_path}/scripts/annotate_gene_mutation_counts.pl \\
  --mutationCounts {mutation_counts} \\
  --vcf {input} > {output} """.format(
            script_path=config.param("DEFAULT", "forge_location"),
            mutation_counts=config.param(ini_section, "mutation_counts_file", type="filepath"),
            input=input,
            output=output
    ))


def homozygosity(input, output, ini_section="homozygosity"):

    return Job(
        [input],
        [output],
        [
            [ini_section, "module_perl"],
            [ini_section, "module_mugqic_tools"]
        ],
        command="""\
{script_path}/scripts/annotate_homozygosity.pl \\
  --vcf {input} > {output}""".format(
            script_path=config.param("DEFAULT", "forge_location"),
            input=input,
            output=output
        ))


def num_private_variants(sample_name, input, output):
    
    return Job(
        [input],
        [output],
        [
            ["get_num_private_variants", "module_perl"],
            ["get_num_private_variants", "module_mugqic_tools"]
        ],
        command="""\
{script_path}/scripts/get_num_private_variants.pl \\
  --vcf {input} \\
  --name {sampleName} > {output}""".format(
            script_path=config.param("DEFAULT", "forge_location"),
            input=input,
            sampleName=sample_name,
            output=output
         ))


def omim(input, output, ini_section="omim"):

    return Job(
        [input],
        [output],
        [
            [ini_section, "module_perl"],
            [ini_section, "module_mugqic_tools"]
        ],
        command="""\
{script_path}/scripts/annotate_omim.pl \\
  --vcf {input} \\
  --omimdb {omimdb} > {output} """.format(
            script_path=config.param("DEFAULT", "forge_location"),
            omimdb=config.param(ini_section, "omimdb", type="filepath"),
            input=input,
            output=output
    ))


def predict_roh(input, output):
    
    return Job(
        [input],
        [output],
        [
            ["predict_roh", "module_perl"],
            ["predict_roh", "module_mugqic_tools"]
        ],
        command="""\
{script_path}/scripts/predict_roh.pl \\
  --vcf {input} \\
  --windowSize 25 \\
  --maxHetsInWindow 2 > {output}""".format(
            script_path=config.param("DEFAULT", "forge_location"),
            input=input,
            output=output
        ))            


def prev_seen(input, output, vardb, prev_seen_details_thresh, ini_section="prev_seen"):
    
    return Job(
        [input, vardb],
        [output],
        [
            [ini_section, "module_perl"],
            [ini_section, "module_mugqic_tools"]
        ],
        command="""\
{script_path}/scripts/annotate_vardb.pl \\
  --variantsDB {vardb} \\
  --vcfFile {input} \\
  --detailsThreshold {prev_seen_details_thresh} > {output}""".format(
            script_path=config.param("DEFAULT", "forge_location"),
            vardb=vardb if vardb else config.param(ini_section, "vardb", type="filepath"),
            input=input,
            prev_seen_details_thresh=prev_seen_details_thresh,
            output=output
        ))


def roh(input, output, roh_file):
    
    return Job(
        [input, roh_file],
        [output, roh_file],
        [
            ["roh", "module_perl"],
            ["roh", "module_mugqic_tools"]
        ],
        command="""\
{script_path}/scripts/annotate_roh.pl \\
  --rohFile {roh_file} \\
  --vcf {input} > {output}""".format(
            script_path=config.param("DEFAULT", "forge_location"),
            roh_file=roh_file,
            input=input,
            output=output
        ))


def ucsc_link(input, output):
    
    return Job(
        [input],
        [output],
        [
            ["ucsc_link", "module_perl"],
            ["ucsc_link", "module_mugqic_tools"]
        ],
        command="""\
{script_path}/scripts/add_hyperlinks.pl \\
  --input {input} > {output}""".format(
            script_path=config.param("DEFAULT", "forge_location"),
            input=input,
            output=output
        ))
            


def update_bam_header(input, output, pgFile, headerFile, ref=None, ini_section="update_bam_header"):

    bwa_ver = config.param("DEFAULT", "module_bwa").split("/")[-1]

    return Job(
        [input, pgFile],
        [headerFile, output, output + ".md5"],
        [
            [ini_section, "module_perl"],
            [ini_section, "module_mugqic_tools"],
	    [ini_section, "module_samtools"],
	    [ini_section, "module_bwa"]
        ],
        command="""\
{script_path}/scripts/updateBAMHeader.pl \\
  --inputFile {input} \\
  --outputFile {output} \\
  --pgFile {pgFile} \\
  --headerFile {headerFile} \\
  --ref_genome {ref} \\
  --bwa_ver {bwa_ver}""".format(
            script_path=config.param("DEFAULT", "forge_location"),
	    input=input,
	    output=output,
            pgFile=pgFile,
	    headerFile=headerFile,
	    ref=ref if ref else config.param(ini_section, "genome_bwa_index", type="filepath"),
	    bwa_ver=bwa_ver
    ))


def update_master_variants(sampleNames, in_vcfs, path_vcf, out_dir, vcf_folder, out_file=None, ini_section="update_master_variants"):

    config_file = config.param(ini_section, "sample_info_file", type="filepath")
    vardb_file = os.path.join(out_dir, config.param(ini_section, "vardb_name"))
    vcflist_file = os.path.join(out_dir, config.param(ini_section, "vcflist_name"))
    other_options = config.param(ini_section, "other_options", required=False)
    in_vcfs.extend([config_file, vardb_file, vcflist_file])

    return Job(
        in_vcfs,
        [vardb_file],
        [
            [ini_section, "module_perl"],
            [ini_section, "module_mugqic_tools"]
        ],
        command="""\
{script_path}/scripts/update_master_variants.pl \\
  --sampleNames {sampleNames} \\
  --sampleConfig {sampleConfig} \\
  --in_vcf_fol {in_vcf_fol} \\
  --vardb {vardb} \\
  --vcflist {vcflist} \\
  --vcffolder {vcffolder}{other_options}{outfile}""".format(
            script_path=config.param("DEFAULT", "forge_location"),
            sampleNames=sampleNames,
            sampleConfig=config_file,
            in_vcf_fol=path_vcf, #os.path.dirname(in_vcfs[0]),
            vardb=vardb_file,
            vcflist=vcflist_file,
            vcffolder=vcf_folder,
            other_options=" \\\n  " + other_options if other_options else "",
            outfile=" > " + out_file if out_file else ""
    ))


def vcf2columns(input, output, cols, escape_cols, col_headers):
    
    return Job(
        [input],
        [output],
        [
            ["vcf2columns", "module_perl"],
            ["vcf2columns", "module_mugqic_tools"]
        ],
        command="""\
{script_path}/scripts/vcf2columns.pl \\
  --vcf {input} \\
  --cols {cols} \\
  --escapeCols {escape_cols} \\
  --headers {col_headers} > {output}""".format(
            script_path=config.param("DEFAULT", "forge_location"),
            input=input,
            cols=cols,
            escape_cols=escape_cols,
            col_headers=col_headers,
            output=output
         ))

