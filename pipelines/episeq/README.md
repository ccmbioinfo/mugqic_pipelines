[TOC]

Episeq Pipeline
===============

The Episeq pipeline takes FASTQ or BAM files (unsorted) as input and produces an differential analysis in the methylome. Currently, only WGSB and RRSB data are supported.

Please ensure that you have Python 2.7.10 available and loaded before executing this pipeline. Otherwise, the program will quit.

In the current working directory, you should have filled out the readset, design and configuration (`.ini`) file. Please refer to the Episeq User Manual for more information and instructions.


Usage
-----
```
#!text

usage: episeq.py [-h] [--help] [-c CONFIG [CONFIG ...]] [-s STEPS]
                 [-o OUTPUT_DIR] [-j {pbs,batch}] [-f] [--report] [--clean]
                 [-l {debug,info,warning,error,critical}] [-d DESIGN]
                 [-r READSETS] [-v]

Version: 2.2.1-beta

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
1- bismark_prepare_genome
2- trim_galore
3- bismark_align
4- picard_merge_sam_files
5- bismark_methylation_caller
6- differential_methylated_pos
7- differential_methylated_regions
```

1- bismark_prepare_genome
-------------------------
This step takes in a reference genome fasta file and convert the sequence for methylation alignment and sequencing. This is a pre-processing step for bismark_align. The step will copy the reference genome to the output directory (if needed), and create the methylome sequence in the directory called `Bisulfite_Genome`, which contains two subdirectories within it. This step only needs to be done once within a project's output folder.
Note that this step can take several hours, depending on the size of the genome.

2- trim_galore
--------------
This step helps improve the alignment efficiency by performing quality trimming via the open source package Trim Galore!. Briefly, this package addresses many of the issues that occur when analyzing RRBS libraries. The pipeline does trimming using only the default options in Trim Galore!. Additional options such as stricter or more relaxed trimming can be entered through the `other_options` parameter in the configuration file. Output files are gzipped `.fq` files. This step can be ignored if the readset has a `.BAM` file.

3- bismark_align
----------------
This step aligns the trimmed reads to a reference genome from `bismark_prepare_genome`. The alignment is done by the open source package Bismark using the default options for RRBS libraries. Additional options can be entered through the "`other_options`" parameter in the configuration file. Output files are `.bam` files, but may be configured to output `cram` or `sam` files. This step can be ignored if the readset contains `.BAM` files

4- picard_merge_sam_files
-------------------------
This step combines all readsets together to form a single `.bam` file for each sample. This is often required when multiple libraries or multiplexing is done before sequencing. Merging at this step allows us to parallelize the pipeline as much as we can before aggregating. `bismark_align` will merge any paired reads into a single file, which takes care some of the work for us.

5- bismark_methylation_caller
-----------------------------
This step extracts the methylation call for every single cytosine analyzed from the Bismark result files.
The following input files are accepted:
    1.	Bismark result files from previous alignment step
    2.	`BAM` files (unsorted) from readset file

6- differential_methylated_pos
------------------------------
This step finds a list of differentially methylated CpG sites with respect to a categorical
phenotype (controls vs. cases). The `BedGraph` files from the previous methylation calling step are first combined
to a `BSRaw` object with the R package `BiSeq`. Then, the `dmpFinder` function from the R package `minfi` is used to
compute a F-test statistic on the beta values for the assayed CpGs in each sample. A p-value is then returned
for each site with the option of correcting them for multiple testing. Differential analysis is done for each
contrast specified in the design file

7- differential_methylated_regions
----------------------------------
This step runs the bumphunting algorithm to locate regions of differential methylation. CpG sites with an average methylation level difference < `delta_beta_threshold` between cases and controls are ignored in the permutation scheme. 
WARNING: This step is slow and requires large amounts of memory!
