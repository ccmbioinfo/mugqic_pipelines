[TOC]


RNAFusion Pipeline
==================
The Gene Fusion pipeline identifies gene fusion events using RNA-seq FASTQ files.  
    
Four separate tools detect fusion events: 
[deFuse](https://sourceforge.net/p/defuse/wiki/DeFuse/), 
[FusionMap](http://www.arrayserver.com/wiki/index.php?title=FusionMap), 
[EricScript](https://sites.google.com/site/bioericscript/home), 
and [INTEGRATE](https://sourceforge.net/p/integrate-fusion/wiki/Home/).

Tophat2 is used to generate precursor files for the INTEGRATE fusion detection tool.

The fusion detection results are combined into one common file (.cff) that gives information about gene fusions including gene names, 
type of fusion (ex. read through vs. gene fusion), and the tools that identified each fusion event. 
Additionally, if DNA sequencing information is available for the samples of interest, 
the Gene Fusion Pipeline can check for DNA support of gene fusions detected from RNA. 

The RNAseq pipeline requires a sampleinfo file to be provided, which is a tab-delimited document with each sample 
as a line with information about sample disease. The first column gives sample name, second column gives the disease name,
third column tells whether the sample comes from a tumor (TP) or normal (NT) tissue. Additional columns can give further 
information about samples.

In addition, a dnabam file must be provided, which gives the name of .bam file(s) associated with RNA-seq sample.
If there is no DNA sequencing associated with the sample, provide the name of an empty file

For validation, an optional file containing fusion gene pairs can be provided with the --valfile flag. Used
only if the validate_fusions step is being run, and assumes that the input sequence data contains the fusions provided in 
the validation file. This step tests the effectiveness of the pipeline in detecting known fusions.

README IS INCOMPLETE BETWEEN THIS POINT...

THIS SUMMARY SECTION IS NOT GENERATED BY THE FUSION PIPELINE. SHOULD A STEP BE ADDED?
Finally, a summary html report is automatically generated by the pipeline at the end of the analysis.
This report contains description
of the sequencing experiment as well as a detailed presentation of the pipeline steps and results.
Various Quality Control (QC) summary statistics are included in the report and additional QC analysis
is accessible for download directly through the report. The report includes also the main references
of the software tools and methods used during the analysis, together with the full list of parameters
that have been passed to the pipeline main script.

An example of the RNA-Seq report for an analysis on Public Corriel CEPH B-cell is available for illustration
purpose only: [RNA-Seq report](http://gqinnovationcenter.com/services/bioinformatics/tools/rnaReport/index.html).

MORE INFORMATION ABOUT THE PIPELINE IS CURRENTLY NOT AVAILABLE. SHOULD THIS BE ADDED?
[Here](https://bitbucket.org/mugqic/mugqic_pipelines/downloads/MUGQIC_Bioinfo_RNA-Seq.pptx) is more
information about the RNA-Seq pipeline that you may find interesting.

... AND THIS POINT.


Usage
-----
```
#!text

usage: rnaseq_fusion.py [-h] [--help] [-c CONFIG [CONFIG ...]] [-s STEPS]
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
1- picard_sam_to_fastq
2- gunzip_fastq
3- merge_fastq
4- defuse
5- fusionmap
6- ericscript
7- tophat2
8- integrate
9- integrate_make_result_file
10- convert_fusion_results_to_cff
11- merge_and_reannotate_cff_fusion
12- check_dna_support_before_next_exon
13- repeat_filter
14- cluster_reann_dnasupp_file
15- fusion_stats
16- validate_fusions
17- delete_fastqs

```
1- picard_sam_to_fastq
----------------------
Convert SAM/BAM files from the input readset file into FASTQ format
if FASTQ files are not already specified in the readset file. Do nothing otherwise.
rerwritten from common.Illumina.picard_sam_to_fastq, make directory for this step under result folder in case the orginal bam file directory is not writtable

2- gunzip_fastq 
---------------
Gunzip .fastq.gz files

3- merge_fastq 
--------------
Merge paired end fastqs of the same sample

4- defuse
---------
Run Defuse to call gene fusions

5- fusionmap 
------------
Run FusionMap to call gene fusions

6- ericscript 
-------------
Run EricScript to call gene fusion

7- tophat2 
----------
Run Tophat2 for Integrate. Determines accepted hits and unmapped reads, and outputs                                                                                            
corresponding .bam files required as input files for integrate step.

8- integrate
------------
Run Integrate to call gene fusions

9- integrate_make_result_file
-----------------------------
Merge infomation from breakpoints.tsv and reads.txt

10- convert_fusion_results_to_cff 
---------------------------------
Convert fusion results of all 4 gene fusion callers to cff format

11- merge_and_reannotate_cff_fusion 
-----------------------------------
Merge all cff files into one single file and reannotate it with given annotation files

12- check_dna_support_before_next_exon 
--------------------------------------
Check DNA support (pair clusters) until the start of next exon/utr

13- repeat_filter 
-----------------
Filter fusions with repetitive boundary sequences by realigning a certain length of sequnces with BWA

14- cluster_reann_dnasupp_file 
------------------------------
Reannotate DNA support (pair clusters) file. This step generates the final category/cluster file,
merged.cff.reann.dnasupp.bwafilter.30.cluster

15- fusion_stats
----------------
Genereates a file containing statistics about detected fusions.

16- validate_fusions
-----------------
Compares the pipeline output in merged.cff.reann.dnasupp.bwafilter.30.cluster with the predetermined fusion gene test file.

17- delete_fastqs 
-----------------
Delete fastqs when all callers' jobs are finished