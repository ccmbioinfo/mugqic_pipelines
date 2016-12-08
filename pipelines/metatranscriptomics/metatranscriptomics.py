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
    * contig_assembly
        Assemble reads into larger contigs,
        map the reads to these contigs
    * search
        Search known microbial sequences using both contigs and singleton reads,
        try bwa, blat, and diamond
    * Further
    Each phase is associated with its own directory
    Within each phase directory, each sample/readset will have its own directory
    """

    # Location for pipeline scripts, to be used by all steps
    # 'metatranscriptomics/scripts'
    script_path = os.path.join(os.path.dirname(__file__), 'scripts')

    def format_fastq_headers(self):
        """
        Mark the headers of the fastq files with /1 or /2 to differentiate the paired-end reads

        Call 'main_add_subID_reads_fastq.pl'
        """
        jobs = []

        output_prefix = 'format_reads'
        for readset in self.readsets:
            output_dir = join(output_prefix, readset.name)

            output1 = join(output_dir, readset.name + '.1.formatted.fastq')
            output2 = join(output_dir, readset.name + '.2.formatted.fastq')

            jobs.append(concat_jobs([mkdir(output1),
                                     mkdir(output2),
                                     Job(input_files=[readset.fastq1, readset.fastq2],
                                         output_files=[output1, output2],
                                         module_entries=[['DEFAULT', 'module_perl']],
                                         command='perl {script_path}/main_add_subID_reads_fastq.pl '
                                                 '{input1} {output1} '
                                                 '{input2} {output2}'.format(script_path=self.script_path,
                                                                             input1=readset.fastq1,
                                                                             output1=output1,
                                                                             input2=readset.fastq2,
                                                                             output2=output2))],
                                    name='format_fastq_headers.' + readset.name))

        return jobs

    def trimmomatic(self):
        jobs = []

        input_prefix = 'format_reads'
        output_prefix = 'format_reads'

        def get_inputs(readset):
            """
            :return: 2 fastq filenames for paired-end reads
            """
            input_dir = join(input_prefix, readset.name)
            return input_dir, \
                   join(input_dir, readset.name + '.1.formatted.fastq'), \
                   join(input_dir, readset.name + '.2.formatted.fastq')

        def get_outputs(readset):
            """
            :return: output directory name,
                     4 fastq filenames
            """
            output_dir = join(output_prefix, readset.name)
            return output_dir, \
                   join(output_dir, readset.name + '.1.qual_paired.fastq'), \
                   join(output_dir, readset.name + '.1.qual_unpaired.fastq'), \
                   join(output_dir, readset.name + '.2.qual_paired.fastq'), \
                   join(output_dir, readset.name + '.2.qual_unpaired.fastq')

        for readset in self.readsets:
            input_dir, input1, input2 = get_inputs(readset)
            output_dir, output_paired1, output_unpaired1, output_paired2, output_unpaired2 = get_outputs(readset)

            job = trimmomatic.trimmomatic(input1,
                                          input2,
                                          output_paired1,
                                          output_unpaired1,
                                          output_paired2,
                                          output_unpaired2,
                                          None,
                                          None,
                                          adapter_file=config.param('trimmomatic', 'adapter_fasta'),
                                          trim_log=join(output_prefix, readset.name, readset.name + '.trim.log'))
            job.name = 'trimmomatic.' + readset.name
            jobs.append(job)

        return jobs

    def merge_overlapping_reads(self):
        """
        Reads from the paired-end fastqs are merged together.

        The merged reads are added to the first paired-end fastq
        """
        jobs = []

        input_prefix = 'format_reads'
        output_prefix = 'format_reads'

        def get_inputs(readset):
            """
            :return: 2 fastq filenames for paired-end reads
            """
            input_dir = join(input_prefix, readset.name)
            return join(input_dir, readset.name + '.1.qual_paired.fastq'), \
                   join(input_dir, readset.name + '.2.qual_paired.fastq')

        def get_flash_params(readset):
            """
            :return: flash's output directory and output prefix
            """
            return join(output_prefix, readset.name), readset.name

        def get_flash_outputs(output_dir, flash_output_prefix):
            """
            Get the filenames that flash will output

            Flash will output 3 files, one for both fastqs, and one for the merged reads
            :return: 2 uncombined fastq files + merged reads fastq file
            """
            return join(output_dir, flash_output_prefix + '.notCombined_1.fastq'), \
                   join(output_dir, flash_output_prefix + '.notCombined_2.fastq'), \
                   join(output_dir, flash_output_prefix + '.extendedFrags.fastq')

        def get_outputs(readset):
            """
            :return: 2 output filenames
            """
            output_dir = join(output_prefix, readset.name)
            return join(output_dir, readset.name + '.1.qual_all.fastq'), \
                   join(output_dir, readset.name + '.2.qual_all.fastq')

        for readset in self.readsets:
            # Get the filenames
            input1, input2 = get_inputs(readset)
            flash_output_dir, flash_output_prefix = get_flash_params(readset)
            flash1, flash2, flash_merged = get_flash_outputs(flash_output_dir, flash_output_prefix)
            output1, output2 = get_outputs(readset)

            # Create jobs
            flash_job = flash.merge_overlapping_reads(input1, input2, flash_output_dir, flash_output_prefix)
            # Put the merged reads into fastq1
            fastq1_job = Job(input_files=[flash1, flash_merged],
                             output_files=[output1],
                             command='cat {flash1} {merged} > {output1}'.format(flash1=flash1, merged=flash_merged,
                                                                                output1=output1))
            # Rename fastq2 to be consistent with fastq1
            fastq2_job = Job(input_files=[flash2],
                             output_files=[output2],
                             command='cp {flash2} {output2}'.format(flash2=flash2, output2=output2))
            jobs.append(concat_jobs([flash_job, fastq1_job, fastq2_job], name='flash.' + readset.name))

        return jobs

    def fastq_to_fasta(self):
        """
        Convert both fastq files to fastas

        Input:
        format_reads/*{1,2}.qual_all.fastq
        Output:
        filter_reads/*{1,2}.qual_all.fasta

        Required since our usearch/5.2.236 only takes fastas as input
        """
        jobs = []

        input_prefix = 'format_reads'
        output_prefix = 'filter_reads'

        for readset in self.readsets:
            input_dir = join(input_prefix, readset.name)
            output_dir = join(output_prefix, readset.name)

            # For both fastq files (paired-end reads)
            for i in (1, 2):
                input_fastq = join(input_dir, readset.name + '.{i}.qual_all.fastq'.format(i=i))
                output_fasta = join(output_dir, readset.name + '.{i}.qual_all.fasta'.format(i=i))
                job_name = 'fastq_to_fasta.{name}.{i}'.format(name=readset.name, i=i)

                jobs.append(seqtk.fastq_to_fasta(input_fastq, output_fasta, name=job_name))

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
                input_fasta = join(input_dir, '{name}.{i}.qual_all.fasta'.format(name=readset.name, i=i))

                output_fasta = join(output_dir, '{name}.{i}.usearch_out.fasta'.format(name=readset.name, i=i))
                output_uc = join(output_dir, '{name}.{i}.usearch_out.uc'.format(name=readset.name, i=i))

                job_name = 'cluster_duplicates.{name}.{i}'.format(name=readset.name, i=i)

                jobs.append(usearch.cluster_duplicates(input_fasta, output_fasta, output_uc, name=job_name))

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
                input_unique_fasta = join(input_unique_dir, '{name}.{i}.usearch_out.fasta'.format(name=readset.name,
                                                                                                  i=i))
                input_uc = join(input_unique_dir, '{name}.{i}.usearch_out.uc'.format(name=readset.name, i=i))

                read_description = join(output_dir, '{name}.{i}.read_description.json'.format(name=readset.name, i=i))
                output_fastq = join(output_dir, '{name}.{i}.unique.fastq'.format(name=readset.name, i=i))
                output_fasta = join(output_dir, '{name}.{i}.unique.fasta'.format(name=readset.name, i=i))

                job_name = 'remove_duplicates.{name}.{i}'.format(name=readset.name, i=i)

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
                jobs.append(infernal.cmscan(rfam_path=config.param('cmscan', 'rfam_location'),
                                            query=join(input_dir,
                                                       '{name}.{i}.unique.fasta'.format(name=readset.name, i=i)),
                                            tblout=join(output_dir,
                                                        '{name}.{i}.infernalout'.format(name=readset.name, i=i)),
                                            log_path=join(output_dir,
                                                          '{name}.{i}.rRNA.log'.format(name=readset.name, i=i)),
                                            name='{step}.{name}'.format(step=self.cmscan.__name__, name=readset.name)))

        return jobs

    def identify_rrna(self):
        """
        Read the output from infernal and determine the FASTQ IDs that are rRNA

        NOTE: this step replaces main_get_infernal_fromfile_1tophit.pl

        Input:
        *.{1,2}.infernalout
        *.{1,2}.read_description.json

        Output:
        *{1,2}.rrna_ids.json            - JSON file w/ rRNA IDs
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

                jobs.append(Job(name='{step}.{readset}.{i}'.format(step=self.identify_rrna.__name__,
                                                                   readset=readset.name,
                                                                   i=i),
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
        *{1,2}.rrna_ids.json        - contains IDs of rRNA reads
        *{1,2}.unique.json          - original reads

        Output:
        *{1,2}.rrna.fastq
        *{1,2}.not_rrna.fastq       - new filtered data set
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
        *.{1,2}.not_rrna.fastq

        Output:
        *.{1,2}.host.sai
        """
        jobs = []

        input_prefix = 'filter_reads'
        output_prefix = 'filter_reads'

        for readset in self.readsets:
            input_dir = join(input_prefix, readset.name)
            output_dir = join(output_prefix, readset.name)

            host_db = config.param(self.align_to_host.__name__, 'host_db')

            # input_fastq[1], input_fastq[2]
            input_fastq = {i: join(input_dir,
                                   '{name}.{i}.not_rrna.fastq'.format(name=readset.name, i=i)) for i in (1, 2)}
            # alignment[1], alignment[2]
            alignment = {i: join(output_dir, '{name}.{i}.host.sai'.format(name=readset.name, i=i)) for i in (1, 2)}

            output_sam = join(output_dir, '{name}.host.sam'.format(name=readset.name))

            for i in (1, 2):
                # bwa aln job
                jobs.append(Job(name='{step}.aln.{readset}.{i}'.format(step=self.align_to_host.__name__,
                                                                       readset=readset.name,
                                                                       i=i),
                                input_files=[input_fastq[i], host_db],
                                output_files=[alignment[i]],
                                module_entries=[[self.align_to_host.__name__, 'module_bwa']],
                                command='bwa aln -t 4 {host_db} {input_fastq}'
                                        '> {alignment}'.format(host_db=host_db,
                                                               input_fastq=input_fastq[i],
                                                               alignment=alignment[i])))

            # bwa sampe job
            jobs.append(Job(name='{step}.sampe.{readset}'.format(step=self.align_to_host.__name__,
                                                                 readset=readset.name),
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
        jobs = []

        input_prefix = 'filter_reads'
        output_prefix = 'filter_reads'

        for readset in self.readsets:
            input_dir = join(input_prefix, readset.name)
            output_dir = join(output_prefix, readset.name)

            input_sam = join(input_dir, '{name}.host.sam'.format(name=readset.name))
            sorted_bam = join(output_dir, '{name}.host_sorted.bam'.format(name=readset.name))
            unmapped_sam = join(output_dir, '{name}.host.bwaout'.format(name=readset.name))
            output_ids = join(output_dir, '{name}.host_ids.json')

            sort_sam_job = Job(name='{step}.sort_host.{name}'.format(step=self.identify_host_reads.__name__,
                                                                     name=readset.name),
                               input_files=[input_sam],
                               output_files=[sorted_bam],
                               module_entries=[[self.identify_host_reads.__name__, 'module_samtools']],
                               command='samtools view -bS {input_sam}'
                                       '| samtools sort -n -o {sorted_bam}'.format(input_sam=input_sam,
                                                                                   sorted_bam=sorted_bam))

            unmapped_reads_job = Job(name='{step}.unmapped_reads.{name}'.format(step=self.identify_host_reads.__name__,
                                                                                name=readset.name),
                                     input_files=[sorted_bam],
                                     output_files=[unmapped_sam],
                                     module_entries=[[self.identify_host_reads.__name__, 'module_sam']],
                                     command='samtools view -F 4 {sorted_bam}'
                                             '> {unmapped_sam}'.format(sorted_bam=sorted_bam,
                                                                       unmapped_sam=unmapped_sam))

            extract_ids_job = Job(name='{step}.extract_IDs.{name}'.format(step=self.identify_host_reads.__name__,
                                                                          name=readset.name),
                                  input_files=[unmapped_sam],
                                  output_files=[output_ids],
                                  module_entries=[[self.identify_host_reads.__name__, 'module_perl']],
                                  command='perl {script_path}/extract_ids_from_sam.py '
                                          '--sam {unmapped_sam} '
                                          '--id-file {output_ids}'.format(script_path=self.script_path,
                                                                          unmapped_sam=unmapped_sam,
                                                                          output_ids=output_ids))

            jobs.extend([sort_sam_job, unmapped_reads_job, extract_ids_job])

        return jobs

    def remove_host_reads(self):
        pass
        # TODO

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
            self.remove_host_reads #12
        ]


if __name__ == '__main__':
    Metatranscriptomics()
