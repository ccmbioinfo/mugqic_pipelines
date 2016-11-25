Episeq Pipeline
===============

The Episeq pipeline takes methylation sequencing data and conducts a differential methylation analysis among the given samples. The pipeline takes FASTQ or BAM files (unsorted) as input and produces an differential analysis in the methylome. Currently, only WGBS and RRSB data are supported.

You should have filled out the readset, design and configuration (`*.ini`) file. Please refer to the Episeq User Manual for more information and instructions.

[TOC]

Quick Start
-----------


Usage
-----
```
#!text

usage: episeq.py [-h] [--help] [-c CONFIG [CONFIG ...]] [-s STEPS]
                 [-o OUTPUT_DIR] [-j {pbs,batch}] [-f] [--report] [--clean]
                 [-l {debug,info,warning,error,critical}] [-d DESIGN]
                 [-r READSETS] [-v]

Version: 2.2.0

For more documentation, visit our website: https://bitbucket.org/mugqic/mugqic_pipelines/

optional arguments:
  -h                    show this help message and exit
  --help                show detailed description of pipeline and steps
  -c CONFIG [CONFIG ...], --config CONFIG [CONFIG ...]
                        config INI-style list of files; config parameters are
                        overwritten based on files order
  -s STEPS, --steps STEPS
                        step range e.g. '1-5', '3,6,7', '2,4-8'
  -o OUTPUT_DIR, --output-dir OUTPUT_DIR
                        output directory (default: current)
  -j {pbs,batch}, --job-scheduler {pbs,batch}
                        job scheduler type (default: pbs)
  -f, --force           force creation of jobs even if up to date (default:
                        false)
  --report              create 'pandoc' command to merge all job markdown
                        report files in the given step range into HTML, if
                        they exist; if --report is set, --job-scheduler,
                        --force, --clean options and job up-to-date status are
                        ignored (default: false)
  --clean               create 'rm' commands for all job removable files in
                        the given step range, if they exist; if --clean is
                        set, --job-scheduler, --force options and job up-to-
                        date status are ignored (default: false)
  -l {debug,info,warning,error,critical}, --log {debug,info,warning,error,critical}
                        log level (default: info)
  -d DESIGN, --design DESIGN
                        design file
  -r READSETS, --readsets READSETS
                        readset file
  -v, --version         show the version information and exit

Steps:
------
 1- bismark_prepare_genome,
 2- pre_qc_check,
 3- trim_galore,
 4- bismark_align,
 5- merge_bismark_alignment_report,
 6- picard_merge_sam_files,
 7- merged_nuc_stats,
 8- bismark_deduplicate,
 9- calc_dedup_nucleotide_coverage,
10- bismark_methylation_caller,
11- bismark_html_report_generator,
12- differential_methylated_pos,
13- differential_methylated_regions
```

## Episeq Pipeline Steps
### 1. bismark_prepare_genome
This step takes in a reference genome fasta file and convert the sequence for methylation alignment and sequencing. This is a pre-processing step for [bismark_align](#4._bismark_align). The step will link the reference genome to the output directory (if needed), and create the methylome sequence in the directory called `Bisulfite_Genome`, which contains two subdirectories within it. Each subdirectory contains a converted sequence of either C->T and G->A conversions. Each of these converted sequences are used to align the bisulfite-treated reads. Since bisulfite sequencing causes unmethylated C to covert to U and later interpreted as T, this step allows alignment to be made without excessive mismatches due to the bisulfite treatment. This step only needs to be done once within a project's output folder.

| Job Attribute | Value |
|:----------|:------|
| Output directory name: | `bismark_prepare_genome` |
| Job name prefix | `bismark_prepare_genome` |
| Requires | None
| Blocks | [bismark_align](#4._bismark_align) <br /> [bismark_methylation_caller](#10._bismark_methylation_caller) |


__Note:__ Depending on the size of the genome, this step can take several hours to complete.

### 2. pre_qc_check
| Job Attribute | Value |
|:----------|:------|
| Output directory name: | `pre_qc_check` |
| Job name prefix | `pre_qc_check` |
| Requires | None |
| Blocks | None |

### 3. trim_galore
This step helps improve the alignment efficiency by performing quality trimming via the open source package Trim Galore!. Briefly, this package addresses many of the issues that occur when analyzing bisulfite sequencing libraries. The pipeline does trimming using only the default options in Trim Galore. Additional options such as stricter or more relaxed trimming can be entered through the `other_options` parameter in the configuration file. Output files are gzipped `.fq` files. This step can be ignored if the readset has a `.BAM` file.

| Job Attribute | Value |
|:----------|:------|
| Output directory name: | `trimmed` |
| Job name prefix | `trim_galore` |
| Requires | None |
| Blocks | [bismark_align](#4._bismark_align) |

### 4. bismark_align
This step aligns the trimmed reads to a reference genome from `bismark_prepare_genome`. The alignment is done by the open source package Bismark using the default stringency settings. By default, the settings can be somewhat strict, but is essential to avoid mismatches from sequencing error. Additional options can be entered through the "`other_options`" parameter in the configuration file. Output files are `.bam` files, but may be configured to output `cram` or `sam` files. This step can be ignored if the readset contains `.BAM` files

| Job Attribute | Value |
|:----------|:------|
| Output directory name: | `aligned` |
| Job name prefix | `bismark_align` |
| Requires | [bismark_prepare_genome](#1._bismark_prepare_genome) <br /> [trim_galore](#3._trim_galore) |
| Blocks | [merge_bismark_alignment_report](#5._merge_bismark_alignment_report) <br /> [picard_merge_sam_files](#6._picard_merge_sam_files) |


### 5. merge_bismark_alignment_report
So far, the pipeline has been handling data at a readset level. However, our analysis requires information on a sample level basis. Thus, we begin to collate all of the information we have collected so for. This step focuses on merging the alignment report that yields information about the mapping quality of the reads.

| Job Attribute | Value |
|:----------|:------|
| Output directory name: | `merged` |
| Job name prefix | `merge_bismark_alignment_report` |
| Requires | [bismark_align](#4._bismark_align) |
| Blocks | [bismark_html_report_generator](#11._bismark_html_report_generator) |

### 6. picard_merge_sam_files
This step combines all readsets together to form a single `.bam` file for each sample. This is often required when multiple libraries or multiplexing is done before sequencing. Merging at this step allows us to parallelize the pipeline as much as we can before aggregating. `bismark_align` will merge any paired reads into a single file, which takes care some of the work for us.

| Job Attribute | Value |
|:----------|:------|
| Output directory name: | `merged` |
| Job name prefix | `picard_merge_sam_files`<br/> `symlink_readset_sample_bam` |
| Requires | [bismark_align](#4._bismark_align) |
| Blocks | [bismark_deduplicate](#8._bismark_deduplicate) <br/> [merged_nuc_stats](#7._merged_nuc_stats) |

### 7. merged_nuc_stats
There isn't a dedicated merge script for the nucleotide coverage information, so we will have to rerun the analysis here to generate a new report file.

| Job Attribute | Value |
|:----------|:------|
| Output directory name: | `merged` |
| Job name prefix | `bam2nuc` |
| Requires | [picard_merge_sam_files](#6._picard_merge_sam_files) |
| Blocks | None |

### 8. bismark_deduplicate
Merging readsets reads can cause all sorts of artifacts and errors. Depending on the experiment, you may have duplicate reads that are frequently found across multiple readsets. Thus, deduplication will help eliminate some background noise that may occur. This yields a processed bam file that is ready for calling and downstream analysis.

| Job Attribute | Value |
|:----------|:------|
| Output directory name: | `dedup` |
| Job name prefix | `bismark_deduplicate`<br/> `skip_rrbs_deduplicate` |
| Requires | [picard_merge_sam_files](#6._picard_merge_sam_files) |
| Blocks | [calc_dedup_nucleotide_coverage](#9._calc_dedup_nucleotide_coverage) <br/> [bismark_methylation_caller](#10._bismark_methylation_caller) <br/> [bismark_html_report_generator](#11._bismark_html_report_generator) |

### 9. calc_dedup_nucleotide_coverage
One more calculation is needed because the bam file has been modified since the last time this calculation has been done. The reason to keep each report file is for diagnostic and troubleshooting purposes. These calculations yield a metric that can be compared as the data moves through the pipeline. For this pipeline, this is the final read coverage information.

| Job Attribute | Value |
|:----------|:------|
| Output directory name: | `dedup` |
| Job name prefix | `bam2nuc` |
| Requires | [bismark_deduplicate](#8._bismark_deduplicate) |
| Blocks | [bismark_html_report_generator](#11._bismark_html_report_generator) |

### 10. bismark_methylation_caller
This step extracts the methylation call for every single cytosine analyzed from the Bismark result files.
The following input files are accepted:
    1.	Bismark result files from previous alignment step
    2.	`BAM` files (unsorted) from readset file

| Job Attribute | Value |
|:----------|:------|
| Output directory name: | `methyl_calls` |
| Job name prefix | `bismark_methylation_caller` |
| Requires | [bismark_prepare_genome](#1._bismark_prepare_genome) <br/> [bismark_deduplicate](#8._bismark_deduplicate) |
| Blocks | [bismark_html_report_generator](#11._bismark_html_report_generator) <br/> [differential_methylated_pos](#12._differential_methylated_pos) <br/>[differential_methylated_regions](13._differential_methylated_regions) |

### 11. bismark_html_report_generator
This job summarizes all data from steps 3-10 into one HTML report file. It contains diagrams that summarizes the various report files Bismark creates in it's proccessing toolkit. This can serve as an excellent overview for the quality of the sample data and could make all other Bismark output reports redundant.

| Job Attribute | Value |
|:----------|:------|
| Output directory name: | `bismark_summary_report` |
| Job name prefix | `bismark_report` |
| Requires | [merge_bismark_alignment_report](#5._merge_bismark_alignment_report) <br/> [bismark_deduplicate](#8._bismark_deduplicate) <br/> [calc_dedup_nucleotide_coverage](#9._calc_dedup_nucleotide_coverage) <br /> [bismark_methylation_caller](#10._bismark_methylation_caller)|
| Blocks | None |

### 12. differential_methylated_pos
This step finds a list of differentially methylated CpG sites with respect to a categorical
phenotype (controls vs. cases). The `BedGraph` files from the previous methylation calling step are first combined
to a `BSRaw` object with the R package `BiSeq`. Then, the `dmpFinder` function from the R package `minfi` is used to
compute a F-test statistic on the beta values for the assayed CpGs in each sample. A p-value is then returned
for each site with the option of correcting them for multiple testing. Differential analysis is done for each
contrast specified in the design file

| Job Attribute | Value |
|:----------|:------|
| Output directory name: | `differential_methylated_pos` |
| Job name prefix | `differential_methylated_pos` |
| Requires | [bismark_methylation_caller](#10._bismark_methylation_caller) |
| Blocks | None |

### 13. differential_methylated_regions
This step runs the bumphunting algorithm to locate regions of differential methylation. CpG sites with an average methylation level difference < `delta_beta_threshold` between cases and controls are ignored in the permutation scheme. 

| Job Attribute | Value |
|:----------|:------|
| Output directory name: | `differential_methylated_regions` |
| Job name prefix | `differential_methylated_regions` |
| Requires | [bismark_methylation_caller](#10._bismark_methylation_caller) |
| Blocks | None |

__WARNING:__ This step is slow and requires large amounts of memory!
