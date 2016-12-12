#!/usr/bin/env python

# Python Standard Modules
import logging
import os.path
from os.path import join
import sys

# Append mugqic_pipelines directory to Python library path
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0])))))

# MUGQIC Modules
from core.config import config
from core.job import Job, concat_jobs, mkdir
from pipelines import common
from bfx import trimmomatic
from bfx import flash
from bfx import seqtk
from bfx import usearch
from bfx import infernal
from bfx import trinity
from bfx import bwa
from bfx import samtools

log = logging.getLogger(__name__)


class Metatranscriptomics(common.Illumina):
    """
    Based off of the pipeline located here:
    http://bioinformatics-ca.github.io/analysis_of_metagenomic_data_module6_lab_2016/

    There are 4 phases to this pipeline so far:
    * format_reads
        Format read headers,
        trim reads,
        merge reads together
    * filter_reads
        Remove duplicates (temporarily),
        remove rRNA,
        remove reads from host
    * contigs
        Assemble reads into larger contigs,
        map the reads to these contigs
    * search
        Search known microbial sequences using both contigs and singleton reads,
        try bwa, blat, and diamond
    * Further phases in progress...

    Each phase is associated with its own directory
    Within each phase directory, each sample/readset will have its own directory
    """

    # Location for pipeline scripts, to be used by all steps
    # 'metatranscriptomics/scripts'
    script_path = os.path.join(os.path.dirname(__file__), 'scripts')

    def format_fastq_headers(self):
        """
        Mark the FASTQ IDs with /1 or /2 to differentiate the paired-end reads

        Call 'main_add_subID_reads_fastq.pl'

        Input:
        <readset fastq1>
        <readset fastq2>

        Output:
        format_reads/*.{1,2}.formatted.fastq
        """
        jobs = []

        output_prefix = 'format_reads'
        for readset in self.readsets:
            output_dir = join(output_prefix, readset.name)

            output_fastq = {}
            for i in (1, 2):
                output_fastq[i] = join(output_dir, '{name}.{i}.formatted.fastq'.format(name=readset.name, i=i))

            jobs.append(concat_jobs([mkdir(output_fastq[1]),
                                     mkdir(output_fastq[2]),
                                     Job(input_files=[readset.fastq1, readset.fastq2],
                                         output_files=[output_fastq[1], output_fastq[2]],
                                         module_entries=[['DEFAULT', 'module_perl']],
                                         command='perl {script_path}/main_add_subID_reads_fastq.pl '
                                                 '{input1} {output1} '
                                                 '{input2} {output2}'.format(script_path=self.script_path,
                                                                             input1=readset.fastq1,
                                                                             output1=output_fastq[1],
                                                                             input2=readset.fastq2,
                                                                             output2=output_fastq[2]))],
                                    name='{step}.{name}'.format(step=self.format_fastq_headers.__name__,
                                                                name=readset.name)))

        return jobs

    def trimmomatic(self):
        """
        Remove adaptors and trim low quality sequences

        Input:
        format_reads/*{1,2}.formatted.fastq

        Output:
        format_reads/*.{1,2}.qual_paired.fastq
        format_reads/*.{1,2}.qual_unpaired.fastq
        """
        jobs = []

        input_prefix = 'format_reads'
        output_prefix = 'format_reads'

        for readset in self.readsets:
            input_dir = join(input_prefix, readset.name)
            output_dir = join(output_prefix, readset.name)

            input_fastq, output_paired, output_unpaired = {}, {}, {}
            for i in (1, 2):
                input_fastq[i] = join(input_dir, '{name}.{i}.formatted.fastq'.format(name=readset.name, i=i))

                output_paired[i] = join(output_dir, '{name}.{i}.qual_paired.fastq'.format(name=readset.name, i=i))
                output_unpaired[i] = join(output_dir, '{name}.{i}.qual_unpaired.fastq'.format(name=readset.name, i=i))

            job = trimmomatic.trimmomatic(input_fastq[1],
                                          input_fastq[2],
                                          output_paired[1],
                                          output_unpaired[1],
                                          output_paired[2],
                                          output_unpaired[2],
                                          None,
                                          None,
                                          adapter_file=config.param(self.trimmomatic.__name__, 'adapter_fasta'),
                                          trim_log=join(output_dir, '{name}.trim.log'.format(name=readset.name)))
            job.name = '{step}.{name}'.format(step=self.trimmomatic.__name__, name=readset.name)
            jobs.append(job)

        return jobs

    def merge_overlapping_reads(self):
        """
        Reads from the paired-end fastqs are merged together using FLASH

        The merged reads are added to the 1st paired-end fastq

        Input:
        format_reads/*.{1,2}.qual_paired.fastq

        Intermediate output:
        format_reads/*.{1,2}.notCombined_1.fastq    - unmerged reads
        format_reads/*.extendedFrags.fastq          - merged reads

        Output:
        format_reads/*.1.qual_all.fastq             - contains fastq1 reads + merged reads
        format_reads/*.2.qual_all.fastq             - contains fastq2 reads
        """
        jobs = []

        input_prefix = 'format_reads'
        output_prefix = 'format_reads'

        for readset in self.readsets:
            input_dir = join(input_prefix, readset.name)
            output_dir = join(output_prefix, readset.name)

            # Get the filenames
            input_fastq, flash_out, output_fastq = {}, {}, {}
            for i in (1, 2):
                input_fastq[i] = join(input_dir, '{name}.{i}.qual_paired.fastq'.format(name=readset.name, i=i))
                flash_out[i] = join(output_dir, '{name}.notCombined_{i}.fastq'.format(name=readset.name, i=i))
                output_fastq[i] = join(output_dir, '{name}.{i}.qual_all.fastq'.format(name=readset.name, i=i))
            flash_merged = join(output_dir, '{name}.extendedFrags.fastq'.format(name=readset.name))

            # Create jobs
            flash_job = flash.merge_overlapping_reads(fastq1=input_fastq[1],
                                                      fastq2=input_fastq[2],
                                                      output_dir=output_dir,
                                                      output_prefix=readset.name)

            # Put the merged reads into fastq1
            fastq1_job = Job(input_files=[flash_out[1], flash_out[2]],
                             output_files=[output_fastq[1]],
                             command='cat {flash1} {merged} > {output1}'.format(flash1=flash_out[1],
                                                                                merged=flash_merged,
                                                                                output1=output_fastq[1]))

            # Rename fastq2 to be consistent with fastq1
            fastq2_job = Job(input_files=[flash_out[2]],
                             output_files=[output_fastq[2]],
                             command='cp {flash2} {output2}'.format(flash2=flash_out[2], output2=output_fastq[2]))

            jobs.append(
                concat_jobs([flash_job,
                             fastq1_job,
                             fastq2_job],
                            name='{step}.{name}'.format(step=self.merge_overlapping_reads.__name__,
                                                        name=readset.name)))

        return jobs

    def fastq_to_fasta(self):
        """
        Convert both fastq files to fastas

        Required since our usearch/5.2.236 only takes fastas as input

        Input:
        format_reads/*{1,2}.qual_all.fastq

        Output:
        filter_reads/*{1,2}.qual_all.fasta
        """
        jobs = []

        input_prefix = 'format_reads'
        output_prefix = 'filter_reads'

        for readset in self.readsets:
            input_dir = join(input_prefix, readset.name)
            output_dir = join(output_prefix, readset.name)

            # For both fastq files (paired-end reads)
            for i in (1, 2):
                jobs.append(
                    seqtk.fastq_to_fasta(
                        fastq=join(input_dir, '{name}.{i}.qual_all.fastq'.format(name=readset.name, i=i)),
                        fasta=join(output_dir, '{name}.{i}.qual_all.fasta'.format(name=readset.name, i=i)),
                        job_name='{step}.{name}'.format(step=self.fastq_to_fasta.__name__, name=readset.name)))

        return jobs

    def cluster_duplicates(self):
        """
        Cluster duplicate reads together

        Input:
        filter_reads/*.{1,2}.qual_all.fasta

        Output:
        filter_reads/*.{1,2}.usearch_out.fasta
                     *.{1,2}.usearch_out.uc         - contains the fasta IDs and cluster IDs
        """
        jobs = []

        input_prefix = 'filter_reads'
        output_prefix = 'filter_reads'

        for readset in self.readsets:
            input_dir = join(input_prefix, readset.name)
            output_dir = join(output_prefix, readset.name)

            # For both paired-end reads
            for i in (1, 2):
                jobs.append(
                    usearch.cluster_duplicates(
                        input_fasta=join(input_dir, '{name}.{i}.qual_all.fasta'.format(name=readset.name, i=i)),
                        output_fasta=join(output_dir, '{name}.{i}.usearch_out.fasta'.format(name=readset.name, i=i)),
                        output_uc=join(output_dir, '{name}.{i}.usearch_out.uc'.format(name=readset.name, i=i)),
                        job_name='{step}.{name}.{i}'.format(step=self.cluster_duplicates.__name__,
                                                            name=readset.name,
                                                            i=i)))

        return jobs

    def remove_duplicates(self):
        """
        Remove duplicate reads but keep track of the number of duplicates for each read

        Input:
        format_reads/*{1,2}.qual_all.fastq          - original duplicated fastq
        filter_reads/*{1,2}.usearch_out.fasta
        filter_reads/*{1,2}.usearch_out.uc

        Output:
        filter_reads/*{1,2}.unique.fastq                - de-duplicated fastq
        filter_reads/*{1,2}.unique.fasta                - de-duplicated fasta
        filter_reads/*{1,2}.read_description.json       - stores ID, number of duplicates, length
        """
        jobs = []

        input_original_prefix = 'format_reads'
        input_unique_prefix = 'filter_reads'
        output_prefix = 'filter_reads'

        for readset in self.readsets:
            input_original_dir = join(input_original_prefix, readset.name)
            input_unique_dir = join(input_unique_prefix, readset.name)
            output_dir = join(output_prefix, readset.name)

            for i in (1, 2):
                input_fastq = join(input_original_dir, '{name}.{i}.qual_all.fastq'.format(name=readset.name, i=i))
                input_unique_fasta = join(input_unique_dir,
                                          '{name}.{i}.usearch_out.fasta'.format(name=readset.name, i=i))
                input_uc = join(input_unique_dir, '{name}.{i}.usearch_out.uc'.format(name=readset.name, i=i))

                read_description = join(output_dir, '{name}.{i}.read_description.json'.format(name=readset.name, i=i))
                output_fastq = join(output_dir, '{name}.{i}.unique.fastq'.format(name=readset.name, i=i))
                output_fasta = join(output_dir, '{name}.{i}.unique.fasta'.format(name=readset.name, i=i))

                job_name = '{step}.{name}.{i}'.format(step=self.remove_duplicates.__name__, name=readset.name, i=i)

                jobs.append(Job(name=job_name,
                                input_files=[input_fastq, input_unique_fasta, input_uc],
                                output_files=[read_description, output_fastq, output_fasta],
                                command='python {script_path}/remove_duplicates.py '
                                        '--input-fastq {input_fastq} '
                                        '--unique-fasta {input_unique_fasta} '
                                        '--unique-uc {input_uc} '
                                        '--output-read-description {read_description} '
                                        '--output-fastq {output_fastq} '
                                        '--output-fasta {output_fasta}'.format(script_path=self.script_path,
                                                                               input_fastq=input_fastq,
                                                                               input_unique_fasta=input_unique_fasta,
                                                                               input_uc=input_uc,
                                                                               read_description=read_description,
                                                                               output_fastq=output_fastq,
                                                                               output_fasta=output_fasta)))

            return jobs

    def cmscan(self):
        """
        Runs cmscan to predict which reads are rRNA

        Note: this is a computationally intensive step

        Input:
        filter_reads/*{1,2}.unique.fasta

        Ouptut:
        filter_reads/*{1,2}.infernalout         - record of which reads are rRNA
        filter_reads/*{1,2].rRNA.log
        """
        jobs = []

        input_prefix = 'filter_reads'
        output_prefix = 'filter_reads'

        for readset in self.readsets:
            input_dir = join(input_prefix, readset.name)
            output_dir = join(output_prefix, readset.name)

            for i in (1, 2):
                jobs.append(
                    infernal.cmscan(
                        rfam_path=config.param('cmscan', 'rfam_location'),
                        query=join(input_dir, '{name}.{i}.unique.fasta'.format(name=readset.name, i=i)),
                        tblout=join(output_dir, '{name}.{i}.infernalout'.format(name=readset.name, i=i)),
                        log_path=join(output_dir, '{name}.{i}.rRNA.log'.format(name=readset.name, i=i)),
                        name='{step}.{name}.{i}'.format(step=self.cmscan.__name__, name=readset.name, i=i)))

        return jobs

    def identify_rrna(self):
        """
        Read the output from infernal and determine the FASTQ IDs that are rRNA

        NOTE: this step replaces main_get_infernal_fromfile_1tophit.pl

        Input:
        filter_reads/*.{1,2}.infernalout
        filter_reds/*.{1,2}.read_description.json

        Output:
        filter_reads/*{1,2}.rrna_ids.json            - JSON file w/ rRNA IDs
        """
        jobs = []

        input_prefix = 'filter_reads'
        output_prefix = 'filter_reads'

        for readset in self.readsets:
            input_dir = join(input_prefix, readset.name)
            output_dir = join(output_prefix, readset.name)

            for i in (1, 2):
                infernalout = join(input_dir, '{name}.{i}.infernalout'.format(name=readset.name, i=i))
                read_description = join(input_dir, '{name}.{i}.read_description.json'.format(name=readset.name, i=i))

                output_ids = join(output_dir, '{name}.{i}.rrna_ids.json'.format(name=readset.name, i=i))

                jobs.append(
                    Job(name='{step}.{readset}.{i}'.format(step=self.identify_rrna.__name__, readset=readset.name, i=i),
                        input_files=[infernalout, read_description],
                        output_files=[output_ids],
                        command='python {script_path}/identify_rrna.py '
                                '--read-description {read_description} '
                                '--infernalout {infernalout} '
                                '--apply-cutoff '
                                '--max-evalue {max_evalue} '
                                '--min-percent-identity {min_percent_identity} '
                                '--out-ids {output_ids}'.format(script_path=self.script_path,
                                                                read_description=read_description,
                                                                infernalout=infernalout,
                                                                max_evalue=0.001,
                                                                min_percent_identity=90,
                                                                output_ids=output_ids)))

        return jobs

    def remove_rrna(self):
        """
        Remove rRNA reads

        Input:
        filter_reads/*{1,2}.rrna_ids.json        - contains IDs of rRNA reads
        filter_reads/*{1,2}.unique.json          - original reads

        Output:
        filter_reads/*{1,2}.rrna.fastq
        filter_reads/*{1,2}.not_rrna.fastq       - new filtered data set
        """
        jobs = []

        input_prefix = 'filter_reads'
        output_prefix = 'filter_reads'

        for readset in self.readsets:
            input_dir = join(input_prefix, readset.name)
            output_dir = join(output_prefix, readset.name)

            for i in (1, 2):
                rrna_ids = join(input_dir, '{name}.{i}.rrna_ids.json'.format(name=readset.name, i=i))
                in_fastq = join(input_dir, '{name}.{i}.unique.fastq'.format(name=readset.name, i=i))

                out_rrna = join(output_dir, '{name}.{i}.rrna.fastq'.format(name=readset.name, i=i))
                out_not_rrna = join(output_dir, '{name}.{i}.not_rrna.fastq'.format(name=readset.name, i=i))

                jobs.append(Job(name='{step}.{readset}'.format(step=self.remove_rrna.__name__, readset=readset.name),
                                input_files=[rrna_ids, in_fastq],
                                output_files=[out_rrna, out_not_rrna],
                                command='python {script_path}/partition_reads_by_id.py '
                                        '--fastq {in_fastq} '
                                        '--id-file {rrna_ids} '
                                        '--included {out_rrna} '
                                        '--excluded {out_not_rrna}'.format(script_path=self.script_path,
                                                                           in_fastq=in_fastq,
                                                                           rrna_ids=rrna_ids,
                                                                           out_rrna=out_rrna,
                                                                           out_not_rrna=out_not_rrna)))

        return jobs

    def align_to_host(self):
        """
        Align reads to a host database

        Input:
        filter_reads/*.{1,2}.not_rrna.fastq

        Output:
        filter_reads/*.{1,2}.host.sai
        """
        jobs = []

        input_prefix = 'filter_reads'
        output_prefix = 'filter_reads'

        for readset in self.readsets:
            input_dir = join(input_prefix, readset.name)
            output_dir = join(output_prefix, readset.name)

            host_db = config.param(self.align_to_host.__name__, 'host_db')

            input_fastq, alignment = {}, {}
            for i in (1, 2):
                input_fastq[i] = join(input_dir, '{name}.{i}.not_rrna.fastq'.format(name=readset.name, i=i))
                alignment[i] = join(output_dir, '{name}.{i}.host.sai'.format(name=readset.name, i=i))

            output_sam = join(output_dir, '{name}.host.sam'.format(name=readset.name))

            for i in (1, 2):
                # bwa aln job
                jobs.append(
                    Job(name='{step}.aln.{name}.{i}'.format(step=self.align_to_host.__name__, name=readset.name, i=i),
                        input_files=[input_fastq[i], host_db],
                        output_files=[alignment[i]],
                        module_entries=[[self.align_to_host.__name__, 'module_bwa']],
                        command='bwa aln -t 4 {host_db} {input_fastq}'
                                '> {alignment}'.format(host_db=host_db,
                                                       input_fastq=input_fastq[i],
                                                       alignment=alignment[i])))

            # bwa sampe job
            jobs.append(
                Job(name='{step}.sampe.{readset}'.format(step=self.align_to_host.__name__, readset=readset.name),
                    input_files=[host_db, alignment[1], alignment[2], input_fastq[1], input_fastq[2]],
                    output_files=[output_sam],
                    module_entries=[[self.align_to_host.__name__, 'module_bwa']],
                    command='bwa sampe {host_db} {alignment1} {alignment2} '
                            '{input_fastq1} {input_fastq2} '
                            '> {merged_sam}'.format(host_db=host_db,
                                                    alignment1=alignment[1],
                                                    alignment2=alignment[2],
                                                    input_fastq1=input_fastq[1],
                                                    input_fastq2=input_fastq[2],
                                                    merged_sam=output_sam)))

        return jobs

    def identify_host_reads(self):
        """
        Identify reads that map to the host's genome

        Input:
        filter_reads/*.host.sam

        Output:
        filter_reads/*.host_ids.json
        """
        jobs = []

        input_prefix = 'filter_reads'
        output_prefix = 'filter_reads'

        for readset in self.readsets:
            input_dir = join(input_prefix, readset.name)
            output_dir = join(output_prefix, readset.name)

            input_sam = join(input_dir, '{name}.host.sam'.format(name=readset.name))
            sorted_bam = join(output_dir, '{name}.host_sorted.bam'.format(name=readset.name))
            host_sam = join(output_dir, '{name}.host.bwaout'.format(name=readset.name))

            output_ids = join(output_dir, '{name}.host_ids.json'.format(name=readset.name))

            # Sort the sam file
            sort_sam_job = Job(
                name='{step}.sort_host.{name}'.format(step=self.identify_host_reads.__name__, name=readset.name),
                input_files=[input_sam],
                output_files=[sorted_bam],
                module_entries=[[self.identify_host_reads.__name__, 'module_samtools']],
                command='samtools view -bS {input_sam}'
                        '| samtools sort -n -o {sorted_bam}'.format(input_sam=input_sam,
                                                                    sorted_bam=sorted_bam))

            # Remove reads that did not map to the host
            remove_unmapped_reads_job = Job(
                name='{step}.remove_unmapped.{name}'.format(step=self.identify_host_reads.__name__, name=readset.name),
                input_files=[sorted_bam],
                output_files=[host_sam],
                module_entries=[[self.identify_host_reads.__name__, 'module_samtools']],
                command='samtools view -F 4 {sorted_bam}'
                        '> {host_sam}'.format(sorted_bam=sorted_bam,
                                              host_sam=host_sam))

            # Given reads that mapped to the host, extract the FASTQ IDs in JSON
            # These are the reads we want to remove
            extract_ids_job = Job(
                name='{step}.extract_IDs.{name}'.format(step=self.identify_host_reads.__name__, name=readset.name),
                input_files=[host_sam],
                output_files=[output_ids],
                command='python {script_path}/extract_ids_from_sam.py '
                        '--sam {unmapped_sam} '
                        '--id-file {output_ids}'.format(script_path=self.script_path,
                                                        unmapped_sam=host_sam,
                                                        output_ids=output_ids))

            jobs.extend([sort_sam_job, remove_unmapped_reads_job, extract_ids_job])

        return jobs

    def remove_host_reads(self):
        """
        Remove the host reads from the data set

        Filter the input FASTQs by using the host IDs in the JSON file

        Input:
        filter_reads/*.host_ids.json         - contains the IDs of host reads
        filter_reads/*.{1,2}.not_rrna.fastq  - original data set including host reads

        Output:
        filter_reads/*.{1,2}.not_host.fastq  - all non-host reads
        filter_reads/*.{1,2}.host.fasq       - all host reads
        """
        jobs = []

        input_prefix = 'filter_reads'
        output_prefix = 'filter_reads'

        for readset in self.readsets:
            input_dir = join(input_prefix, readset.name)
            output_dir = join(output_prefix, readset.name)

            for i in (1, 2):
                id_file = join(input_dir, '{name}.host_ids.json'.format(name=readset.name))
                input_fastq = join(input_dir, '{name}.{i}.not_rrna.fastq'.format(name=readset.name, i=i))

                output_not_host = join(output_dir, '{name}.{i}.not_host.fastq'.format(name=readset.name, i=i))
                output_host = join(output_dir, '{name}.{i}.host.fastq'.format(name=readset.name, i=i))

                jobs.append(
                    Job(name='{step}.{name}.{i}'.format(step=self.remove_host_reads.__name__, name=readset.name, i=i),
                        input_files=[id_file, input_fastq],
                        output_files=[output_host, output_not_host],
                        command='python {script_path}/partition_reads_by_id.py '
                                '--fastq {input_fastq} '
                                '--id-file {id_file} '
                                '--included {output_host} '
                                '--excluded {output_not_host}'.format(script_path=self.script_path,
                                                                      input_fastq=input_fastq,
                                                                      id_file=id_file,
                                                                      output_host=output_host,
                                                                      output_not_host=output_not_host)))

        return jobs

    def return_duplicates(self):
        """
        Return the duplicate reads back to the data set

        Assists with assembly coverage

        Input:
        filter_reads/*.{1,2}.read_description.json
        filter_reads/*.{1,2}.not_host.fastq

        Output:
        filter_reads/*.{1,2}.mRNA.fastq
        """
        jobs = []

        input_prefix = 'filter_reads'
        output_prefix = 'filter_reads'

        for readset in self.readsets:
            for i in (1, 2):
                input_dir = join(input_prefix, readset.name)
                output_dir = join(output_prefix, readset.name)

                input_fastq = join(input_dir, '{name}.{i}.not_host.fastq'.format(name=readset.name, i=i))
                read_description = join(input_dir, '{name}.{i}.read_description.json'.format(name=readset.name, i=i))

                output_fastq = join(output_dir, '{name}.{i}.mRNA.fastq'.format(name=readset.name, i=i))

                jobs.append(
                    Job(name='{step}.{name}.{i}'.format(step=self.return_duplicates.__name__, name=readset.name, i=i),
                        input_files=[input_fastq, read_description],
                        output_files=[output_fastq],
                        command='python {script_path}/return_duplicates.py '
                                '--input-fastq {input_fastq} '
                                '--read-description {read_description} '
                                '--output-fastq {output_fastq}'.format(script_path=self.script_path,
                                                                       input_fastq=input_fastq,
                                                                       read_description=read_description,
                                                                       output_fastq=output_fastq)))

        return jobs

    def trinity(self):
        """
        Assemble the reads into contigs for searching

        Input:
        filter_reads/*{1,2}.mRNA.fastq

        Output:
        contigs/*.contigs.fasta
        """
        jobs = []

        input_prefix = 'filter_reads'
        output_prefix = 'contigs'

        for readset in self.readsets:
            input_dir = join(input_prefix, readset.name)
            output_dir = join(output_prefix, readset.name)
            trinity_out_dir = join(output_dir, 'trinity_out_dir')

            input_mrna_1 = join(input_dir, '{name}.1.mRNA.fastq'.format(name=readset.name))
            input_mrna_2 = join(input_dir, '{name}.2.mRNA.fastq'.format(name=readset.name))
            input_str = '--left {mrna_1} --right {mrna_2}'.format(mrna_1=input_mrna_1, mrna_2=input_mrna_2)

            trinity_fasta = join(trinity_out_dir, 'Trinity.fasta')
            output_fasta = join(output_dir, '{name}.contigs.fasta'.format(name=readset.name))

            # Create 'Trinity.fasta'
            trinity_job = concat_jobs([mkdir(trinity_fasta),
                                       trinity.trinity(input_files=[input_mrna_1, input_mrna_2],
                                                       trinity_fasta=trinity_fasta,
                                                       output_directory=trinity_out_dir,
                                                       reads_option=input_str)])

            # Rename output to '*.contigs.fasta'
            rename_job = Job(input_files=[trinity_fasta],
                             output_files=[output_fasta],
                             command='cp {trinity_fasta} {output_fasta}'.format(trinity_fasta=trinity_fasta,
                                                                                output_fasta=output_fasta))

            jobs.append(concat_jobs([trinity_job, rename_job],
                                    name='{step}.{name}'.format(step=self.trinity.__name__, name=readset.name)))

        return jobs

    def index_contigs(self):
        """
        Index the assembled contigs with 'bwa index' and 'samtools faidx'

        Input:
        contigs/*.contigs.fasta

        Output:
        contigs/*.contigs.fasta.amb     - bwa index
        contigs/*.contigs.fasta.sa      - bwa index
        contigs/*.contigs.fasta.ann     - bwa index
        contigs/*.contigs.fasta.pac     - bwa index
        contigs/*.contigs.fasta.bwt     - bwa index
        contigs/*.contigs.fasta.fai     - samtools faidx
        """
        jobs = []

        contig_prefix = 'contigs'

        for readset in self.readsets:
            contig_dir = join(contig_prefix, readset.name)

            contigs = join(contig_dir, '{name}.contigs.fasta'.format(name=readset.name))

            # 'bwa index'
            bwa_index_job = bwa.index(contigs)
            bwa_index_job.name = '{step}.bwa_index.{name}'.format(step=self.index_contigs.__name__, name=readset.name)

            # 'samtools faidx'
            samtools_faidx_job = samtools.faidx(contigs)
            samtools_faidx_job.name = '{step}.samtools_faidx.{name}'.format(step=self.index_contigs.__name__,
                                                                            name=readset.name)
            jobs.extend([bwa_index_job, samtools_faidx_job])

        return jobs

    @property
    def steps(self):
        return [
            self.format_fastq_headers,
            self.trimmomatic,
            self.merge_overlapping_reads,  # 3
            self.fastq_to_fasta,
            self.cluster_duplicates,
            self.remove_duplicates,  # 6
            self.cmscan,
            self.identify_rrna,
            self.remove_rrna,  # 9
            self.align_to_host,
            self.identify_host_reads,
            self.remove_host_reads,  # 12
            self.return_duplicates,
            self.trinity,
            self.index_contigs,  # 15
        ]


if __name__ == '__main__':
    Metatranscriptomics()
