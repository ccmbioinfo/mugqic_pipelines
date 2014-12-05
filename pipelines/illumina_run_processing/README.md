# Illumina Run Processing Pipeline


Version: 2.0.0-beta


## Overview


The standard MUGQIC Illumina Run Processing pipeline uses the Illumina bcl2fastq
software to convert and demultiplex the base call files to fastq files. The
pipeline runs some QCs on the raw data, on the fastq and on the alignment.

## Sample Sheets
The pipeline uses two input sample sheets. The first one is the standard Casava
sheet, a csv file having the following columns (please refer to the Illumina
Casava user guide):

- `SampleID`
- `FCID`
- `SampleRef`
- `Index`
- `Description`
- `Control`
- `Recipe`
- `Operator`
- `SampleProject`

The second sample sheet is called the Nanuq run sheet. It's a csv file with the
following minimal set of mandatory columns (the column order in the file doesn't
matter)

- `ProcessingSheetId` Must be the same as the `SampleID` from the Casava Sheet.
- `Name` The sample name put in RG headers of bam files and on filename on disk.
- `Run` The run number.
- `Region` The lane number.
- `Library Barcode` The library barcode put in .bam's RG headers and on disk
- `Library Source` The type of library. If this value contains 'RNA' or 'cDNA',
STAR will be used to make the aligmnent, otherwise, bwa\_mem will be used
- `BED Files` The name of the BED file containing the genomic targets. This is
the 'filename' parameter passed to the 'fetch\_bed\_file\_command'


On this page:

[TOC]


## Usage
```
#!bash
mugqic_pipeline/pipelines/illumina_run_processing/illumina_run_processing.py --help
```


## Steps

 This pipeline performs the following steps:


### Step 1: index


Generate a file with all the indexes found in the index-reads of the run.

The input barcode file is a two columns tsv file. Each line has a
"barcode\_sequence" and the corresponding "barcode\_name". This file can be
generated by a LIMS.

The output is a tsv file named "RUNFOLDER\_LANENUMBER.metrics" that will be
saved in the output directory. This file has four columns, the barcode/index
sequence, the index name, the number of reads and the number of reads that have
passed the filter.


### Step 2: fastq


Launch fastq generation from Illumina raw data using BCL2FASTQ conversion
software.

The index base mask is calculated according to the sample and run configuration;
and also according the mask parameters received (first/last index bases). The
Casava sample sheet is generated with this mask. The default number of
mismatches allowed in the index sequence is 1 and can be overrided with an
command line argument. No demultiplexing occurs when there is only one sample in
the lane.

An optional notification command can be launched to notify the start of the
fastq generation with the calculated mask.


### Step 3: align


Align the reads from the fastq file, sort the resulting .bam and create an index
of that .bam.

An basic aligment is performed on a sample when the "SampleRef" field of the
Illumina sample sheet match one of the regexp in the configuration file and the
corresponding genome (and indexes) are installed.

STAR is used as a splice-junctions aware aligner when the sample
"library\_source" is cDNA or contains "RNA"; otherwise BWA\_mem is used to align 
the reads.


### Step 4: picard\_mark\_duplicates


Runs Picard mark duplicates on the sorted bam file.


### Step 5: metrics


This step runs a series of multiple metrics collection jobs and the output bam
from mark duplicates.

- Picard CollectMultipleMetrics: A collection of picard metrics that runs at the
same time to save on I/O.
    - CollectAlignmentSummaryMetrics,
    - CollectInsertSizeMetrics,
    - QualityScoreDistribution,
    - MeanQualityByCycle,
    - CollectBaseDistributionByCycle
- BVATools DepthOfCoverage: Using the specified "BED Files" in the sample sheet,
calculate the coverage of each target region.
- Picard CalculateHsMetrics: Calculates a set of Hybrid Selection specific
metrics from the BAM file. The bait and interval list is automatically created
from the specicied "BED Files".


### Step 6: blast


Run blast on a subsample of the reads of each sample to find the 20 most
frequent hits.

The "runBlast.sh" util from MUGQIC Tools is used. The number of reads to
subsample can be configured by sample or for the whole lane. The output will be
in the "Blast\_sample" folder, under the Unaligned folder.


### Step 7: qc\_graphs


Generate some QC Graphics and a summary XML file for each sample using 
[BVATools](https://bitbucket.org/mugqic/bvatools/).

Files are created in a 'qc' subfolder of the fastq directory. Examples of
output graphic:

- Per cycle qualities, sequence content and sequence length;
- Known sequences (adaptors);
- Abundant Duplicates;


### Step 8: md5


Create md5 checksum files for the fastq, bam and bai using the system 'md5sum'
util.

One checksum file is created for each file.


### Step 9: start\_copy\_notification


Send an optional notification for the processing completion.

The command used is in the configuration file. This step is skipped when no
command is provided.


### Step 10: copy


Copy processed files to another place where they can be served or loaded into a
LIMS.

The destination folder and the command used can be set in the configuration
file.


### Step 11: end\_copy\_notification


Send an optional notification to notify that the copy is finished.

The command used is in the configuration file. This step is skipped when no
command is provided.

