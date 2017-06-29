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


#---------------------------------------------STEP1 STARTS HERE----------------------------------------------
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

    # TODO: use common/trimmomatic
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
                                          adapter_file=config.param(self.trimmomatic.__name__, 'adapter_fasta',
                                                                    type='filepath'),
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
#---------------------------------------------STEP1 ENDS HERE----------------------------------------------

#---------------------------------------------STEP2 STARTS HERE----------------------------------------------
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
        filter_reads/*{1,2}.read_description.txt       - stores ID, number of duplicates, length
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

                read_description = join(output_dir, '{name}.{i}.read_description.txt'.format(name=readset.name, i=i))
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
#---------------------------------------------STEP2 ENDS HERE----------------------------------------------


#---------------------------------------------STEP3 STARTS HERE----------------------------------------------
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
                        database=config.param('cmscan', 'rfam_location', type='filepath'),
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
        filter_reds/*.{1,2}.read_description.txt

        Output:
        filter_reads/*{1,2}.rrna_ids.txt            - txt file w/ rRNA IDs
        """
        jobs = []

        input_prefix = 'filter_reads'
        output_prefix = 'filter_reads'

        for readset in self.readsets:
            input_dir = join(input_prefix, readset.name)
            output_dir = join(output_prefix, readset.name)

            for i in (1, 2):
                infernalout = join(input_dir, '{name}.{i}.infernalout'.format(name=readset.name, i=i))
                read_description = join(input_dir, '{name}.{i}.read_description.txt'.format(name=readset.name, i=i))

                output_ids = join(output_dir, '{name}.{i}.rrna_ids.txt'.format(name=readset.name, i=i))

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
        filter_reads/*{1,2}.rrna_ids.txt        - contains IDs of rRNA reads
        filter_reads/*{1,2}.unique.fastq          - original reads

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
                rrna_ids = join(input_dir, '{name}.{i}.rrna_ids.txt'.format(name=readset.name, i=i))
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
                                        '--excluded {out_not_rrna} '
                                        '--out-format fastq'.format(script_path=self.script_path,
                                                                    in_fastq=in_fastq,
                                                                    rrna_ids=rrna_ids,
                                                                    out_rrna=out_rrna,
                                                                    out_not_rrna=out_not_rrna)))

        return jobs
#---------------------------------------------STEP3 ENDS HERE----------------------------------------------


#---------------------------------------------STEP4 STARTS HERE----------------------------------------------
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

            host_db = config.param(self.align_to_host.__name__, 'host_db', type='filepath')

            input_fastq, alignment = {}, {}
            for i in (1, 2):
                input_fastq[i] = join(input_dir, '{name}.{i}.not_rrna.fastq'.format(name=readset.name, i=i))
                alignment[i] = join(output_dir, '{name}.{i}.host.sai'.format(name=readset.name, i=i))

                # bwa aln job
                jobs.append(
                    bwa.aln(
                        query=host_db,
                        target=input_fastq[i],
                        output=alignment[i],
                        num_threads=4,
                        name='{step}.aln.{name}.{i}'.format(step=self.align_to_host.__name__, name=readset.name, i=i)))

            output_sam = join(output_dir, '{name}.host.sam'.format(name=readset.name))

            # bwa sampe job
            jobs.append(
                bwa.sampe(target=host_db,
                          sai1=alignment[1],
                          sai2=alignment[2],
                          fastq1=input_fastq[1],
                          fastq2=input_fastq[2],
                          output=output_sam,
                          name='{step}.sampe.{name}'.format(step=self.align_to_host.__name__, name=readset.name)))

        return jobs

    def identify_host_reads(self):
        """
        Identify reads that map to the host's genome

        Input:
        filter_reads/*.host.sam

        Output:
        filter_reads/*.host_ids.txt
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

            output_ids = join(output_dir, '{name}.host_ids.txt'.format(name=readset.name))

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
        filter_reads/*.host_ids.txt         - contains the IDs of host reads
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
                id_file = join(input_dir, '{name}.host_ids.txt'.format(name=readset.name))
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
                                '--excluded {output_not_host} '
                                '--out-format fastq'.format(script_path=self.script_path,
                                                            input_fastq=input_fastq,
                                                                id_file=id_file,
                                                                output_host=output_host,
                                                                output_not_host=output_not_host)))

        return jobs
#---------------------------------------------STEP4 ENDS HERE----------------------------------------------

#---------------------------------------------STEP5 STARTS HERE----------------------------------------------
    def return_duplicates(self):
        """
        Return the duplicate reads back to the data set

        Assists with assembly coverage

        Input:
        filter_reads/*.{1,2}.read_description.txt
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
                read_description = join(input_dir,
                        '{name}.{i}.read_description.txt'.format(name=readset.name, i=i))

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
#---------------------------------------------STEP5 ENDS HERE----------------------------------------------

#---------------------------------------------STEP6 STARTS HERE----------------------------------------------
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
        contigs/*.contigs.fasta.bwt     - bwa index
        contigs/*.contigs.fasta.fai     - samtools faidx
        """
        jobs = []

        contig_prefix = 'contigs'

        for readset in self.readsets:
            contig_dir = join(contig_prefix, readset.name)

            contigs = join(contig_dir, '{name}.contigs.fasta'.format(name=readset.name))

            # 'bwa index'
            bwa_index_job = bwa.index(contigs, options=config.param('bwa_mem', 'bwa_index_options'))
            bwa_index_job.name = '{step}.bwa_index.{name}'.format(step=self.index_contigs.__name__, name=readset.name)

            # 'samtools faidx'
            samtools_faidx_job = samtools.faidx(contigs)
            samtools_faidx_job.name = '{step}.samtools_faidx.{name}'.format(step=self.index_contigs.__name__,
                                                                            name=readset.name)
            jobs.extend([bwa_index_job, samtools_faidx_job])

        return jobs

    def align_to_contigs(self):
        """
        Align both sets of reads to the assembled contigs

        Input:
        filter_reads/*.{1,2}.mRNA.fastq
        contigs/*.contigs.fasta

        Output:
        contigs/*.trinity.sam
        """
        jobs = []

        input_reads_prefix = 'filter_reads'
        input_contigs_prefix = 'contigs'

        output_prefix = 'contigs'

        for readset in self.readsets:
            input_reads_dir = join(input_reads_prefix, readset.name)
            input_contigs_dir = join(input_contigs_prefix, readset.name)
            output_dir = join(output_prefix, readset.name)

            contigs = join(input_contigs_dir, '{name}.contigs.fasta'.format(name=readset.name))

            input_reads, sai = {}, {}
            for i in (1, 2):
                input_reads[i] = join(input_reads_dir, '{name}.{i}.mRNA.fastq'.format(name=readset.name, i=i))
                sai[i] = join(output_dir, '{name}.{i}.trinity.sai'.format(name=readset.name, i=i))

                # bwa aln job
                jobs.append(
                    bwa.aln(query=contigs,
                            target=input_reads[i],
                            output=sai[i],
                            num_threads=4,
                            name='{step}.aln.{name}.{i}'.format(step=self.align_to_contigs.__name__,
                                                                name=readset.name,
                                                                i=i)))

            output_sam = join(output_dir, '{name}.trinity.sam'.format(name=readset.name))

            # bwa sampe job
            jobs.append(
                bwa.sampe(target=contigs,
                          sai1=sai[1],
                          sai2=sai[2],
                          fastq1=input_reads[1],
                          fastq2=input_reads[2],
                          output=output_sam,
                          name='{step}.sampe.{name}'.format(step=self.align_to_contigs.__name__, name=readset.name)))

        return jobs

    def identify_contigs_reads(self):
        """
        Identify reads that map to the contigs

        Input: contigs/cow_readset.trinity.sam

        Output: contigs/cow_readset.trinity_id.txt
        """
        jobs = []

        input_prefix = 'contigs'
        output_prefix = 'contigs'

        for readset in self.readsets:
            input_dir = join(input_prefix, readset.name)
            output_dir = join(output_prefix, readset.name)

            input_sam = join(input_dir,'{name}.trinity.sam'.format(name=readset.name))
            sorted_bam = join(output_dir,'{name}.trinity.bam'.format(name=readset.name))
            contigs_sam = join(output_dir,'{name}.trinity.bwaout'.format(name=readset.name))
            id_file = join(output_dir,'{name}.trinity_id.txt'.format(name=readset.name))

            # sort the sam file
            sort_sam_job = Job(
                    name='{step}.sort_contigs.{name}'.format(step=self.identify_contigs_reads.__name__,
                        name=readset.name),
                    input_files=[input_sam],
                    output_files=[sorted_bam],
                    module_entries=[[self.identify_contigs_reads.__name__,'module_samtools']],
                    command='samtools view -bS {input_sam} | samtools sort -n'
                    ' -o {sorted_bam}'.format(input_sam=input_sam,
                        sorted_bam=sorted_bam))

            # remove reads that did not map to the contigs
            remove_unmapped_reads_job= Job(
                    name='{step}.remove_unmapped.{name}'.format(step=self.identify_contigs_reads.__name__,
                        name=readset.name),
                    input_files=[sorted_bam],
                    output_files=[contigs_sam],
                    module_entries=[[self.identify_contigs_reads.__name__,'module_samtools']],
                    command='samtools view -F 4 {sorted_bam} > '
                    '{contigs_sam}'.format(sorted_bam=sorted_bam,
                        contigs_sam=contigs_sam))

            # Extract IDs from .bwaout file in JSON format
            extract_ID_job = Job(
                    name='{step}.extract_ids.{name}'.format(step=self.extract_singletons.__name__,
                        name=readset.name),
                    input_files=[contigs_sam],
                    output_files=[id_file],
                    command='python {script_path}/extract_ids_from_sam.py '
                        '--sam {input_sam} '
                        '--id-file {id_file}'.format(script_path=self.script_path,
                                                input_sam=contigs_sam,
                                                id_file=id_file))

            jobs.extend([sort_sam_job, remove_unmapped_reads_job, extract_ID_job])

        return jobs

    def extract_singletons(self):
        """
        Extract singletons and mapped putative mRNA

        Input:
        {readset}.{i}.mRNA.fastq

        Output:
        {readset}.{i}.singleton.fastq
        {readset}.{i}.mappedreads.fastq
        """
        jobs = []

        input_filter_prefix = 'filter_reads'
        input_contigs_prefix = 'contigs'
        output_prefix = 'contigs'

        for readset in self.readsets:
            input_filter_dir = join(input_filter_prefix, readset.name)
            input_contigs_dir = join(input_contigs_prefix, readset.name)
            output_dir = join(output_prefix, readset.name)

            id_file = join(input_contigs_dir,'{name}.trinity_id.txt'.format(name=readset.name))

            output_singleton = {}
            for i in (1,2):
                # Extract singleton reads from fastq
                input_fastq = join(input_filter_dir,
                        '{name}.{i}.mRNA.fastq'.format(name=readset.name, i=i))
                output_singleton[i] = join(output_dir,
                        '{name}.{i}.singletons.fastq'.format(name=readset.name, i=i))
                output_mRNA_mappedreads = join(output_dir,
                        '{name}.{i}.mRNA_mappedreads.fastq'.format(name=readset.name, i=i))

                jobs.append(
                        Job(name='{step}.{name}.{i}.get_singletons'.format(step=self.extract_singletons.__name__,
                                                                    name=readset.name, i=i),
                            input_files=[id_file, input_fastq],
                            output_files=[output_singleton[i],output_mRNA_mappedreads],
                            command='python {script_path}/partition_reads_by_id.py '
                                    '--fastq {input_fastq} '
                                    '--id-file {id_file} '
                                    '--included {output_mRNA_mappedreads} '
                                    '--excluded {output_singleton} '
                                    '--out-format fastq'.format(script_path=self.script_path,
                                                                input_fastq=input_fastq,
                                                                id_file=id_file,
                                                                output_mRNA_mappedreads=output_mRNA_mappedreads,
                                                                output_singleton=output_singleton[i])))

            # Merge two singleton files into one, and extract length file
            merged = join(output_dir,
                    '{name}.singletons.fastq'.format(name=readset.name))
            length_file = join(output_dir,'{name}.singletons_length.txt'.format(name=readset.name))

            jobs.append(Job(
                name='{step}.{name}.merge_fastq'.format(step=self.extract_singletons.__name__,
                                                        name=readset.name),
                input_files=[output_singleton[1], output_singleton[2]],
                output_files=[merged],
                command='cat {fastq1} {fastq2} > {merged}'.format(fastq1=output_singleton[1],
                                                                fastq2=output_singleton[2],
                                                                merged=merged)))

            jobs.append(
                    Job(name='{step}.{name}.get_length'.format(step=self.extract_singletons.__name__,
                                                                    name=readset.name),
                    input_files=[merged],
                    output_files=[length_file],
                    command='python {script_path}/main_get_sequence_length.py '
                            '--fastq {fastq} '
                            '--output {length_file}'.format(script_path=self.script_path,
                                                            fastq=merged,
                                                            length_file=length_file)))
        return jobs

    def get_mapping_table(self):
        """
        Generates txt file to find the number of reads that are mapped to
        contigs

        Input:
        {readset}.trinity_id.txt
        {readset}.contigs.fasta

        Output:
        {readset}.contigs.IDs_length.txt
        """
        jobs = []

        input_prefix = 'contigs'
        output_prefix = 'contigs'

        for readset in self.readsets:
            input_dir = join(input_prefix, readset.name)
            output_dir = join(output_prefix, readset.name)

            # Generate maptable to show the contigs and the number of reads
            # thar are mapped to the contigs
            pair_id = join(input_dir, '{name}.trinity_id.txt'.format(name=readset.name))
            input_fasta = join(input_dir,'{name}.contigs.fasta'.format(name=readset.name))
            maptable = join(output_dir,
                    '{name}.contigs.IDs_length.txt'.format(name=readset.name))

            jobs.append(Job(
                    name='{step}.{name}'.format(step=self.get_mapping_table.__name__,
                                                name=readset.name),
                    input_files=[pair_id],
                    output_files=[maptable],
                    command='python {script_path}/get_maptable.py '
                            '--id-file {pair_id} '
                            '--contig-fasta {fasta} '
                            '--maptable {maptable}'.format(script_path=self.script_path,
                                                                pair_id=pair_id,
                                                                fasta=input_fasta,
                                                                maptable=maptable)))

        return jobs
#---------------------------------------------STEP6 ENDS HERE-----------------------------------------

#---------------------------------------------STEP7 STARTS HERE-----------------------------------------
    def bwa_align_contigs(self):
        """
        Align contigs.fasta (single-ended read) using BWA

        Input:
        {readset}.contigs.fasta

        Output:
        {readset}.bwa.sam
        """
        jobs = []

        input_prefix = "contigs"
        output_prefix = "contigs"

        for readset in self.readsets:
            input_dir = join(input_prefix, readset.name)
            output_dir = join(output_prefix, readset.name)

            blast_db = config.param(self.bwa_align_contigs.__name__, 'blast_db', type='filepath')

            input_fasta = join(input_dir, '{name}.contigs.fasta'.format(name=readset.name))
            output_sai = join(output_dir, '{name}.contigs.sai'.format(name=readset.name))

            # bwa aln; alignment to get sai file
            jobs.append(
                    bwa.aln(
                        query=blast_db,
                        target=input_fasta,
                        output=output_sai,
                        num_threads=4,
                        name='{step}.aln.{name}'.format(step=self.bwa_align_contigs.__name__,
                                                        name=readset.name)))

            # bwa samse; converting sai to sam (for single-ended file)
            output_sam = join(output_dir, '{name}.contigs.sam'.format(name=readset.name))

            jobs.append(
                    bwa.samse(
                        target=blast_db,
                        sai=output_sai,
                        fasta=input_fasta,
                        output=output_sam,
                        name='{step}.samse.{name}'.format(step=self.bwa_align_contigs.__name__,
                                                            name=readset.name)))
        return jobs

    def bwa_identify_contigs(self):
        """
        Map microbial genes to the identified contigs using BWA

        Input:
        {readset}.contigs.sam

        Output:
        {readset}.contig.bam
        {readset}.contigs.micro_cds.bwaout
        """
        jobs = []

        input_prefix = 'contigs'
        output_prefix = 'contigs'

        for readset in self.readsets:
            input_dir = join(input_prefix, readset.name)
            output_dir = join(output_prefix, readset.name)

            input_sam = join(input_dir, '{name}.contigs.sam'.format(name=readset.name))
            output_bam = join(output_dir, '{name}.contigs.bam'.format(name=readset.name))
            mapped_sam = join(output_dir, '{name}.contigs.micro_cds.bwaout'.format(name=readset.name))

            # Convert sam to binary bam file, and sort it
            sort_sam_job = Job(
                    name='{step}.sort_sam.{name}'.format(step=self.bwa_identify_contigs.__name__,
                                                        name=readset.name),
                    input_files=[input_sam],
                    output_files=[output_bam],
                    module_entries=[[self.bwa_identify_contigs.__name__, 'module_samtools']],
                    command='samtools view -bS {input_sam} | '
                            'samtools sort -n -o {output_bam}'.format(input_sam=input_sam,
                                                                    output_bam=output_bam))

            # Map microbial genes 
            map_sam_job = Job(
                    name='{step}.map_sam.{name}'.format(step=self.bwa_identify_contigs.__name__,
                                                        name=readset.name),
                    input_files=[output_bam],
                    output_files=[mapped_sam],
                    module_entries=[[self.bwa_identify_contigs.__name__, 'module_samtools']],
                    command='samtools view -F 4 {output_bam} > {mapped_sam}'.format(
                                                                                output_bam=output_bam,
                                                                                mapped_sam=mapped_sam))
            jobs.extend([sort_sam_job, map_sam_job])
        return jobs

    def bwa_contigs_select_reads(self):
        """
        Extract microbial and non-microbial genes from BWA search output

        Input:
        {readset}.contigs.micro_cds.bwaout
        {readset}.contigs.fasta

        Output:
        {readset}.contigs.micro_cds_id_bwa.txt
        {readset}.contigs.micro_cds.fasta
        {readset}.contigs.n_micro_cds.fasta
        """
        jobs = []

        input_prefix = 'contigs'
        output_prefix = 'contigs'

        for readset in self.readsets:
            input_dir = join(input_prefix, readset.name)
            output_dir = join(output_prefix, readset.name)

            input_sam = join(input_dir,'{name}.contigs.micro_cds.bwaout'.format(name=readset.name))
            output_id = join(output_dir,'{name}.contigs.micro_cds_id_bwa.txt'.format(name=readset.name))

            # Extract microbial contigs IDs 
            jobs.append(Job(
                name='{step}.{name}.micro_cds_id'.format(step=self.bwa_contigs_select_reads.__name__,
                                                        name=readset.name),
                input_files=[input_sam],
                output_files=[output_id],
                command='python {script_path}/extract_ids_from_sam.py '
                        '--sam {input_sam} '
                        '--id-file {output_id}'.format(script_path=self.script_path,
                                                        input_sam=input_sam,
                                                        output_id=output_id)))

            # Extract microbial and non-microbial genes in fasta 
            input_fasta = join(input_dir,'{name}.contigs.fasta'.format(name=readset.name))
            output_micro_cds = join(output_dir,'{name}.contigs.micro_cds.fasta'.format(name=readset.name))
            output_n_micro_cds = join(output_dir,'{name}.contigs.n_micro_cds.fasta'.format(name=readset.name))

            jobs.append(Job(
                name='{step}.{name}.select_reads'.format(step=self.bwa_contigs_select_reads.__name__,
                                                        name=readset.name),
                input_files=[input_fasta, output_id],
                output_files=[output_micro_cds, output_n_micro_cds],
                command='python {script_path}/partition_reads_by_id.py '
                        '--fasta {input_fasta} '
                        '--id-file {output_id} '
                        '--included {output_micro_cds} '
                        '--excluded {output_n_micro_cds} '
                        '--out-format fasta'.format(script_path=self.script_path,
                                                    input_fasta=input_fasta,
                                                    output_id=output_id,
                                                    output_micro_cds=output_micro_cds,
                                                    output_n_micro_cds=output_n_micro_cds)))
        return jobs

    def bwa_align_singletons(self):
        """
        Align putative mRNA using BWA 

        Input:
        {readset}.{1,2}.singletons.fastq

        Output:
        {readset}.{1,2}.singletons.sai
        {readset}.singletons.sam
        """
        jobs = []

        input_prefix = 'contigs'
        output_prefix = 'contigs'

        for readset in self.readsets:
            input_dir = join(input_prefix, readset.name)
            output_dir = join(output_prefix, readset.name)
            blast_db = config.param(self.bwa_align_singletons.__name__, 'blast_db', type='filepath')

            input_fastq, alignment = {}, {}
            for i in (1,2):
                # Align putative mRNA 
                input_fastq[i] = join(input_dir,'{name}.{i}.singletons.fastq'.format(name=readset.name,i=i))
                alignment[i] = join(output_dir,'{name}.{i}.singletons.sai'.format(name=readset.name,i=i))

                jobs.append(
                        bwa.aln(
                            query=blast_db,
                            target=input_fastq[i],
                            output=alignment[i],
                            num_threads=4,
                            name='{step}.aln.{name}.{i}'.format(step=self.bwa_align_singletons.__name__,
                                                                name=readset.name,
                                                                i=i)))

            output_sam = join(output_dir,'{name}.singletons.sam'.format(name=readset.name))
            # Generate sam file 
            jobs.append(
                    bwa.sampe(
                        target=blast_db,
                        sai1=alignment[1],
                        sai2=alignment[2],
                        fastq1=input_fastq[1],
                        fastq2=input_fastq[2],
                        output=output_sam,
                        name='{step}.sampe.{name}'.format(step=self.bwa_align_singletons.__name__,
                                                        name=readset.name)))
        return jobs

    def bwa_identify_singletons(self):
        """
        Map microbial database to identified singletons using BWA 

        Input:
        {readset}.singletons.sam

        Output:
        {readset}.singletons.micro_cds.bam
        {readset}.singletons.micro_cds.bwaout
        """
        jobs = []

        input_prefix = 'contigs'
        output_prefix = 'contigs'

        for readset in self.readsets:
            input_dir = join(input_prefix, readset.name)
            output_dir = join(output_prefix, readset.name)

            input_sam = join(input_dir,'{name}.singletons.sam'.format(name=readset.name))
            sorted_bam = join(output_dir,'{name}.singletons.micro_cds.bam'.format(name=readset.name))

            # Sort Sam and convert to binary format
            jobs.append(Job(
                name='{step}.sort_sam.{name}'.format(step=self.bwa_identify_singletons.__name__,
                                                    name=readset.name),
                input_files=[input_sam],
                output_files=[sorted_bam],
                module_entries=[[self.bwa_identify_singletons.__name__,'module_samtools']],
                command='samtools view -bS {input_sam} | '
                        'samtools sort -n -o {sorted_bam}'.format(input_sam=input_sam,
                                                                sorted_bam=sorted_bam)))
            # Map microbial genes to identified singletons
            mapped_sam = join(output_dir,'{name}.singletons.micro_cds.bwaout'.format(name=readset.name))

            jobs.append(Job(
                name='{step}.map_sam.{name}'.format(step=self.bwa_identify_singletons.__name__,
                                                    name=readset.name),
                input_files=[sorted_bam],
                output_files=[mapped_sam],
                module_entries=[[self.bwa_identify_singletons.__name__,'module_samtools']],
                command='samtools view -F 4 {sorted_bam} > {mapped_sam}'.format(sorted_bam=sorted_bam,
                                                                                mapped_sam=mapped_sam)))
        return jobs

    def bwa_singletons_select_reads(self):
        """
        Select microbial and non-microbial gene and generate fasta for each

        Input:
        {readset}.singletons.micro_cds.bwaout
        {readset}.{i}.singletons.fastq

        Output:
        {readset}.singletons.micro_cds_id_bwa.txt
        {readset}.{1,2}.singletons.micro_cds.fasta
        {readset}.{1,2}.singletons.n_micro_cds.fasta
        """
        jobs = []

        input_prefix = 'contigs'
        output_prefix = 'contigs'

        for readset in self.readsets:
            input_dir = join(input_prefix, readset.name)
            output_dir = join(output_prefix, readset.name)

            # Extract IDs from Sam output
            input_sam = join(input_dir,'{name}.singletons.micro_cds.bwaout'.format(name=readset.name))
            output_id = join(output_dir,'{name}.singletons.micro_cds_id_bwa.txt'.format(name=readset.name))

            jobs.append(Job(
                name='{step}.{name}.micro_cds_id'.format(step=self.bwa_singletons_select_reads.__name__,
                                                        name=readset.name),
                input_files=[input_sam],
                output_files=[output_id],
                command='python {script_path}/extract_ids_from_sam.py '
                        '--sam {input_sam} '
                        '--id-file {output_id}'.format(script_path=self.script_path,
                                                        input_sam=input_sam,
                                                        output_id=output_id)))
            # Select microbial and non-microbial genes from singletons
            for i in (1,2):
                input_fastq = join(input_dir,'{name}.{i}.singletons.fastq'.format(name=readset.name,i=i))
                output_micro_cds = join(output_dir,'{name}.{i}.singletons.micro_cds.fasta'.format(name=readset.name,i=i))
                output_n_micro_cds = join(output_dir,'{name}.{i}.singletons.n_micro_cds.fasta'.format(name=readset.name,i=i))

                jobs.append(Job(
                    name='{step}.{i}.{name}.select_read'.format(step=self.bwa_singletons_select_reads.__name__,
                                                                name=readset.name,
                                                                i=i),
                    input_files=[input_fastq, output_id],
                    output_files=[output_micro_cds, output_n_micro_cds],
                    command='python {script_path}/partition_reads_by_id.py '
                            '--fastq {input_fastq} '
                            '--id-file {output_id} '
                            '--included {output_micro_cds} '
                            '--excluded {output_n_micro_cds} '
                            '--out-format fasta'.format(script_path=self.script_path,
                                                        input_fastq=input_fastq,
                                                        output_id=output_id,
                                                        output_micro_cds=output_micro_cds,
                                                        output_n_micro_cds=output_n_micro_cds)))

        return jobs

    def blat_search_contigs(self):
        """
        Search contigs against microbial database using BLAT

        Input:
        {readset}.contigs.n_micro_cds.fasta

        Output:
        {readset}.contigs_1.blatout
        {readset}.contigs_2.blatout
        {readset}.contigs.n_micro_cds.blatout
        """
        jobs = []

        input_prefix = 'contigs'
        output_prefix = 'contigs'

        for readset in self.readsets:
            input_dir = join(input_prefix, readset.name)
            output_dir = join(output_prefix, readset.name)

            blast_db_1 = config.param(self.blat_search_contigs.__name__,'blast_db_1', type='filepath')
            blast_db_2 = config.param(self.blat_search_contigs.__name__,'blast_db_2', type='filepath')

            # Align contigs against microbial genome database using BLAT
            # Due to large size of the database, run BLAT twice with divided
            # databases to avoid out of memory issue
            input_fasta = join(input_dir,'{name}.contigs.n_micro_cds.fasta'.format(name=readset.name))

            output_fasta1 = join(output_dir,'{name}.contigs_1.blatout'.format(name=readset.name))
            output_fasta2 = join(output_dir,'{name}.contigs_2.blatout'.format(name=readset.name))

            jobs.append(Job(
                name='{step}.{name}.blat_1'.format(step=self.blat_search_contigs.__name__,
                                                    name=readset.name),
                input_files=[input_fasta, blast_db_1],
                output_files=[output_fasta1],
                module_entries=[[self.blat_search_contigs.__name__, 'module_blat']],
                command='blat -noHead -minIdentity={minIdentity} '
                        '-minScore={minScore} {blast_db_1} {input_fasta} '
                        '-fine -q={query_type} -t={db_type} '
                        '-out=blast8 {output_fasta1}'.format(minIdentity=90,
                                                            minScore=50,
                                                            blast_db_1=blast_db_1,
                                                            input_fasta=input_fasta,
                                                            query_type='rna',
                                                            db_type='dna',
                                                            output_fasta1=output_fasta1)))

            jobs.append(Job(
                name='{step}.{name}.blat_2'.format(step=self.blat_search_contigs.__name__,
                                                    name=readset.name),
                input_files=[input_fasta, blast_db_2],
                output_files=[output_fasta2],
                module_entries=[[self.blat_search_contigs.__name__, 'module_blat']],
                command='blat -noHead -minIdentity={minIdentity} '
                        '-minScore={minScore} {blast_db_2} {input_fasta} '
                        '-fine -q={query_type} -t={db_type} '
                        '-out=blast8 {output_fasta2}'.format(minIdentity=90,
                                                            minScore=50,
                                                            blast_db_2=blast_db_2,
                                                            input_fasta=input_fasta,
                                                            query_type='rna',
                                                            db_type='dna',
                                                            output_fasta2=output_fasta2)))

            # Merge two alignment files into one
            output_merged = join(output_dir,'{name}.contigs.n_micro_cds.blatout'.format(name=readset.name))

            jobs.append(Job(
                name='{step}.{name}.merge_files'.format(step=self.blat_search_contigs.__name__,
                                                        name=readset.name),
                input_files=[output_fasta1,output_fasta2],
                output_files=[output_merged],
                command='cat {output_fasta1} {output_fasta2} > '
                        '{output_merged}'.format(output_fasta1=output_fasta1,
                                                output_fasta2=output_fasta2,
                                                output_merged=output_merged)))
        return jobs

    def blat_search_singletons(self):
        """
        Search singletons against microbial database using BLAT

        Input:
        {readset}.{1,2}.singletons.n_micro_cds.fasta

        Output:
        {readset}.{1,2}.singletons_1.blatout
        {readset}.{1,2}.singletons_2.blatout
        {readset}.{1,2}.singletons.n_micro_cds.blatout
        {readset}.singletons.n_micro_cds.blatout
        """
        jobs = []

        input_prefix = 'contigs'
        output_prefix = 'contigs'

        for readset in self.readsets:
            input_dir = join(input_prefix, readset.name)
            output_dir = join(output_prefix, readset.name)

            # Align singletons against microbial genome database
            blast_db_1 = config.param(self.blat_search_singletons.__name__,'blast_db_1', type='filepath')
            blast_db_2 = config.param(self.blat_search_singletons.__name__,'blast_db_2', type='filepath')

            output_merged = {}
            for i in (1,2):
                input_fasta = join(input_dir,'{name}.{i}.singletons.n_micro_cds.fasta'.format(name=readset.name,i=i))
                output_fasta1 = join(output_dir,'{name}.{i}.singletons_1.blatout'.format(name=readset.name,i=i))
                output_fasta2 = join(output_dir,'{name}.{i}.singletons_2.blatout'.format(name=readset.name,i=i))

                jobs.append(Job(
                    name='{step}.{name}.{i}.blat_1'.format(step=self.blat_search_singletons.__name__,
                                                            name=readset.name,
                                                            i=i),
                    input_files=[input_fasta, blast_db_1],
                    output_files=[output_fasta1],
                    module_entries=[[self.blat_search_singletons.__name__, 'module_blat']],
                    command='blat -noHead -minIdentity={minIdentity} '
                            '-minScore={minScore} {blast_db_1} {input_fasta} '
                            '-fine -q={query_type} -t={db_type} '
                            '-out=blast8 {output_fasta1}'.format(minIdentity=90,
                                                                minScore=50,
                                                                blast_db_1=blast_db_1,
                                                                input_fasta=input_fasta,
                                                                query_type='rna',
                                                                db_type='dna',
                                                                output_fasta1=output_fasta1)))

                jobs.append(Job(
                    name='{step}.{name}.{i}.blat_2'.format(step=self.blat_search_singletons.__name__,
                                                            name=readset.name,
                                                            i=i),
                    input_files=[input_fasta, blast_db_2],
                    output_files=[output_fasta2],
                    module_entries=[[self.blat_search_singletons.__name__, 'module_blat']],
                    command='blat -noHead -minIdentity={minIdentity} '
                            '-minScore={minScore} {blast_db_2} {input_fasta} '
                            '-fine -q={query_type} -t={db_type} '
                            '-out=blast8 {output_fasta2}'.format(minIdentity=90,
                                                                minScore=50,
                                                                blast_db_2=blast_db_2,
                                                                input_fasta=input_fasta,
                                                                query_type='rna',
                                                                db_type='dna',
                                                                output_fasta2=output_fasta2)))

                # Merge two BLAT results into one
                output_merged[i] = join(output_dir,'{name}.{i}.singletons.n_micro_cds.blatout'.format(name=readset.name,i=i))

                jobs.append(Job(
                    name='{step}.{name}.{i}.merge_files'.format(step=self.blat_search_singletons.__name__,
                                                                name=readset.name,
                                                                i=i),
                    input_files=[output_fasta1, output_fasta2],
                    output_files=[output_merged[i]],
                    command='cat {output_fasta1} {output_fasta2} > '
                            '{output_merged}'.format(output_fasta1=output_fasta1,
                                                    output_fasta2=output_fasta2,
                                                    output_merged=output_merged[i])))

            # Combine two singleton files into one
            output_combined = join(output_dir,'{name}.singletons.n_micro_cds.blatout'.format(name=readset.name))

            jobs.append(Job(
                name='{step}.{name}.combine'.format(step=self.blat_search_singletons.__name__,
                                                    name=readset.name),
                input_files=[output_merged[1], output_merged[2]],
                output_files=[output_combined],
                command='cat {output1} {output2} > {output_combined}'.format(output1=output_merged[1],
                                                                            output2=output_merged[2],
                                                                            output_combined=output_combined)))
        return jobs

    def process_contigs(self):
        """
        Process BLAT output; sort, extract tophits, and generate fasta

        Input:
        {readset}.contigs.n_micro_cds.blatout
        {readset}.contigs.IDs_length.txt
        {readset}.contigs.n_micro_cds.fasta

        Output:
        {readset}.contigs.n_micro_cds_sorted.blatout
        {readset}.contigs.n_micro_cds_pairs.txt
        {readset}.contigs.n_micro_cds_id_blat.txt
        {readset}.contigs.n_micro_cds_blat.fasta
        {readset}.contigs.n_micro_cds_rest.fasta
        """
        jobs = []

        input_prefix = 'contigs'
        output_prefix = 'contigs'

        for readset in self.readsets:
            input_dir = join(input_prefix, readset.name)
            output_dir = join(output_prefix, readset.name)

            blatout = join(input_dir,
                    '{name}.contigs.n_micro_cds.blatout'.format(name=readset.name))

            sort_blat = join(output_dir,
                    '{name}.contigs.n_micro_cds_sorted.blatout'.format(name=readset.name))

            # Sort blat file
            jobs.append(Job(
                name='{step}.{name}.sort'.format(step=self.process_contigs.__name__,
                                                name=readset.name),
                input_files=[blatout],
                output_files=[sort_blat],
                command='python {script_path}/sort_blastout.py '
                        '--infile {blatout} '
                        '--sort {sort_blat} '
                        '--cutoff 10'.format(script_path=self.script_path,
                                            blatout=blatout,
                                            sort_blat=sort_blat)))

            # Process sorted blat file
            processed_blat = join(output_dir,'{name}.contigs.n_micro_cds_pairs.txt'.format(name=readset.name))
            length_file = join(input_dir,'{name}.contigs.IDs_length.txt'.format(name=readset.name))
            id_file = join(output_dir,
                    '{name}.contigs.n_micro_cds_id_blat.txt'.format(name=readset.name))

            jobs.append(Job(
                name='{step}.{name}.get_blast'.format(step=self.process_contigs.__name__,
                                                    name=readset.name),
                input_files=[length_file, sort_blat],
                output_files=[processed_blat, id_file],
                command='python {script_path}/get_blast_1tophit.py '
                        '--length {length_file} '
                        '--infile {sort_blat} '
                        '--outfile {processed_blat} '
                        '--id-file {id_file} '
                        '--cutoff-type 1 '
                        '--cutoff0 100 '
                        '--cutoff1 85 '
                        '--cutoff2 65 '
                        '--cutoff3 60'.format(script_path=self.script_path,
                                            length_file=length_file,
                                            sort_blat=sort_blat,
                                            processed_blat=processed_blat,
                                            id_file=id_file)))

            # Select reads from processed blat file
            input_fasta = join(input_dir,
                    '{name}.contigs.n_micro_cds.fasta'.format(name=readset.name))
            blat_fasta = join(output_dir,
                    '{name}.contigs.n_micro_cds_blat.fasta'.format(name=readset.name))
            rest_fasta = join(output_dir,
                    '{name}.contigs.n_micro_cds_rest.fasta'.format(name=readset.name))

            jobs.append(Job(
                name='{step}.{name}.select_reads'.format(step=self.process_contigs.__name__,
                                                        name=readset.name),
                input_files=[input_fasta, id_file],
                output_files=[blat_fasta, rest_fasta],
                command='python {script_path}/partition_reads_by_id.py '
                        '--fasta {input_fasta} '
                        '--id-file {id_file} '
                        '--included {blat_fasta} '
                        '--excluded {rest_fasta} '
                        '--out-format fasta'.format(script_path=self.script_path,
                                                    input_fasta=input_fasta,
                                                    id_file=id_file,
                                                    blat_fasta=blat_fasta,
                                                    rest_fasta=rest_fasta)))
        return jobs

    def process_singletons(self):
        """
        Process BLAT output; sort, extract tophits, and generate fasta

        Input:
        {readset}.1.singletons.n_micro_cds.blatout
        {readset}.2.singletons.n_micro_cds.blatout
        {readset}.{1,2}.singletons.n_micro_cds.fasta

        Output:
        {readset}.singletons.n_micro_cds_sorted.blatout
        {readset}.singletons.n_micro_cds_pairs.txt
        {readset}.singletons.n_micro_cds_IDs.txt
        {readset}.{1,2}.singletons.n_micro_cds_blat.fasta
        {readset}.{1,2}.singletons.n_micro_cds_rest.fasta
        """
        jobs = []

        input_prefix = 'contigs'
        output_prefix = 'contigs'

        for readset in self.readsets:
            input_dir = join(input_prefix, readset.name)
            output_dir = join(output_prefix, readset.name)

            # merge individual blatout files and sort the combined file
            singleton_blat_1 = join(input_dir,
                    '{name}.1.singletons.n_micro_cds.blatout'.format(name=readset.name))
            singleton_blat_2 = join(input_dir,
                    '{name}.2.singletons.n_micro_cds.blatout'.format(name=readset.name))
            combined_blat = join(output_dir,
                    '{name}.singletons.n_micro_cds.blatout'.format(name=readset.name))

            sorted_blat = join(output_dir,
                    '{name}.singletons.n_micro_cds_sorted.blatout'.format(name=readset.name))

            jobs.append(concat_jobs([
                Job(input_files=[singleton_blat_1, singleton_blat_2],
                    output_files=[combined_blat],
                    command='cat {blat1} {blat2} > {combined}'.format(blat1=singleton_blat_1,
                                                                    blat2=singleton_blat_2,
                                                                    combined=combined_blat)),
                Job(input_files=[combined_blat],
                    output_files=[sorted_blat],
                    command='python {script_path}/sort_blastout.py '
                            '--infile {combined} '
                            '--sort {sort} '
                            '--cutoff 10'.format(script_path=self.script_path,
                                                combined=combined_blat,
                                                sort=sorted_blat))],
                name='{step}.{name}.sort_blat'.format(step=self.process_singletons.__name__,
                                                    name=readset.name)))

            # get top hits from the sorted blast file
            combined_length = join(output_dir,
                    '{name}.singletons_length.txt'.format(name=readset.name))

            tophits_blat = join(output_dir,
                    '{name}.singletons.n_micro_cds_pairs.txt'.format(name=readset.name))
            id_file = join(output_dir,
                '{name}.singletons.n_micro_cds_id_blat.txt'.format(name=readset.name))

            jobs.append(Job(
                    name='{step}.{name}.get_tophits'.format(step=self.process_singletons.__name__,
                                                            name=readset.name),
                    input_files=[combined_length, sorted_blat],
                    output_files=[tophits_blat, id_file],
                    command='python {script_path}/get_blast_1tophit.py '
                            '--length {length} '
                            '--infile {sort} '
                            '--outfile {tophits} '
                            '--id-file {id_file} '
                            '--cutoff-type 0 '
                            '--cutoff0 100 '
                            '--cutoff1 85 '
                            '--cutoff2 65 '
                            '--cutoff3 60'.format(script_path=self.script_path,
                                                length=combined_length,
                                                sort=sorted_blat,
                                                tophits=tophits_blat,
                                                id_file=id_file)))

            # Select reads 
            for i in (1,2):
                input_fasta = join(input_dir,
                        '{name}.{i}.singletons.n_micro_cds.fasta'.format(name=readset.name,
                                                                        i=i))
                blat_fasta = join(output_dir,
                        '{name}.{i}.singletons.n_micro_cds_blat.fasta'.format(name=readset.name,
                                                                            i=i))
                rest_fasta = join(output_dir,
                        '{name}.{i}.singletons.n_micro_cds_rest.fasta'.format(name=readset.name,
                                                                            i=i))
                jobs.append(Job(
                        name='{step}.{i}.{name}.select_reads'.format(step=self.process_singletons.__name__,
                                                                    name=readset.name,
                                                                    i=i),
                        input_files=[input_fasta, id_file],
                        output_files=[blat_fasta, rest_fasta],
                        command='python {script_path}/partition_reads_by_id.py '
                                '--fasta {input_fasta} '
                                '--id-file {id_file} '
                                '--included {blat} '
                                '--excluded {rest} '
                                '--out-format fasta'.format(script_path=self.script_path,
                                                            input_fasta=input_fasta,
                                                            id_file=id_file,
                                                            blat=blat_fasta,
                                                            rest=rest_fasta)))
        return jobs

    def diamond_align_contigs(self):
        """
        Align contigs against non-redundant protein database using DIAMOND

        Input:
        {readset}.contigs.n_micro_cds_rest.fasta

        Output:
        {readset}.contigs.nr.matches.daa
        {readset}.contigs.nr.diamondout

        """
        jobs = []

        input_prefix = 'contigs'
        output_prefix = 'contigs'

        for readset in self.readsets:
            input_dir = join(input_prefix, readset.name)
            output_dir = join(output_prefix, readset.name)
            tmp_dir = join(output_dir, 'diamond_tmp')

            blast_db = config.param(self.diamond_align_contigs.__name__,
                    'blast_db', type='filepath')
            tmp = join(tmp_dir, 'dir_created')
            # Dummy file; Diamond modifies tmp directory as it is running
            # This makes the latest modification time of the directory more
            # recent than generated inputs. As a result, isUp2Date function 
            # in Job.py complains even though all the outputs are up to date
            # To prevent this, create a dummy file
            dummy = join(output_dir, 'dummy')

            # Let database to be the input requirement for making dummy file
            jobs.append(Job(name='{step}.{name}.make_dummy'.format(step=self.diamond_align_contigs.__name__,
                                                                name=readset.name),
                            input_files=[blast_db],
                            output_files=[dummy],
                            command='touch {dummy}'.format(dummy=dummy)))


            input_fasta = join(input_dir,
                    '{name}.contigs.n_micro_cds_rest.fasta'.format(name=readset.name))
            aligned = join(output_dir,
                    '{name}.contigs.nr.matches'.format(name=readset.name))
            # Dummy output file format. Diamond adds .daa to its output
            # If "aligned" is used, mugqic pipeline complains because it can't
            # find its file name. Same for singletons case
            aligned1 = join(output_dir,
                    '{name}.contigs.nr.matches.daa'.format(name=readset.name))

            tmp_dir_job = Job(name='{step}.{name}.Make_tmp_dir'.format(step=self.diamond_align_contigs.__name__,
                                                                        name=readset.name),
                                input_files=[dummy],
                                output_files=[tmp_dir],
                                command='mkdir -p {tmp}'.format(tmp=tmp_dir))

            tmp_file_job = Job(name='{step}.{name}.make_tmp_file'.format(step=self.diamond_align_contigs.__name__,
                                                                        name=readset.name),
                                input_files=[dummy],
                                output_files=[tmp],
                                command='touch {tmp}'.format(tmp=tmp))

            jobs.extend([tmp_dir_job, tmp_file_job])

            jobs.append(Job(
                name='{step}.{name}.align'.format(step=self.diamond_align_contigs.__name__,
                                                    name=readset.name),
                input_files=[tmp, blast_db, input_fasta],
                output_files=[aligned1],
                module_entries=[[self.diamond_align_contigs.__name__,
                    'module_diamond']],
                command='diamond blastx -p 8 -d {blast_db} -q {input_fasta} '
                        '-a {aligned} -t {tmp_dir} -e 10 -k 10'.format(
                                                                    blast_db=blast_db,
                                                                    input_fasta=input_fasta,
                                                                    aligned=aligned,
                                                                    tmp_dir=tmp_dir)))

            # Map read using diamond view
            mapped = join(output_dir,
                    '{name}.contigs.nr.diamondout'.format(name=readset.name))

            jobs.append(Job(
                name='{step}.{name}.map_reads'.format(step=self.diamond_align_contigs.__name__,
                                                        name=readset.name),
                input_files=[aligned1],
                output_files=[mapped],
                module_entries=[[self.diamond_align_contigs.__name__,
                    'module_diamond']],
                command='diamond view -a {aligned} -o {mapped} -f tab'.format(
                                                                            aligned=aligned,
                                                                            mapped=mapped)))
        return jobs

    def diamond_align_singletons(self):
        """
        Align singletons against non-redundant protein database using DIAMOND

        Input:
        {readset}.{1,2}.singletons.n_micro_cds_rest.fasta

        Output:
        {readset}.{1,2}.singletons.nr.matches.daa
        {readset}.singletons.nr.diamondout
        """
        jobs = []

        input_prefix = 'contigs'
        output_prefix = 'contigs'

        for readset in self.readsets:
            input_dir = join(input_prefix, readset.name)
            output_dir = join(output_prefix, readset.name)
            diamond_dir = join(output_dir, 'diamond_tmp')
            check = join(diamond_dir, 'dir_created')

            blast_db = config.param(self.diamond_align_singletons.__name__,
                    'blast_db', type='filepath')

            mapped = {}
            for i in (1,2):
                input_fasta = join(input_dir,
                        '{name}.{i}.singletons.n_micro_cds_rest.fasta'.format(name=readset.name,i=i))
                aligned = join(output_dir,
                        '{name}.{i}.singletons.nr.matches'.format(name=readset.name,i=i))
                aligned1 = join(output_dir,
                        '{name}.{i}.singletons.nr.matches.daa'.format(name=readset.name,i=i))

                jobs.append(Job(
                    name='{step}.{i}{name}.align'.format(step=self.diamond_align_singletons.__name__,
                                                        name=readset.name,i=i),
                    input_files=[check, input_fasta, blast_db],
                    output_files=[aligned1],
                    module_entries=[[self.diamond_align_singletons.__name__,
                        'module_diamond']],
                    command='diamond blastx -p 8 -d {blast_db} -q {input_fasta} '
                            '-a {aligned} -t {tmp_dir} -e 10 -k 10'.format(
                                                                        blast_db=blast_db,
                                                                        input_fasta=input_fasta,
                                                                        aligned=aligned,
                                                                        tmp_dir=diamond_dir)))
                # Map reads using diamond view
                mapped[i] = join(output_dir,
                        '{name}.{i}.singletons.nr.diamondout'.format(name=readset.name,i=i))

                jobs.append(Job(
                    name='{step}.{i}.{name}.map_reads'.format(step=self.diamond_align_singletons.__name__,
                                                            name=readset.name,i=i),
                    input_files=[aligned1],
                    output_files=[mapped[i]],
                    module_entries=[[self.diamond_align_singletons.__name__,
                        'module_diamond']],
                    command='diamond view -a {aligned} -o {mapped} -f tab'.format(
                                                                                aligned=aligned,
                                                                                mapped=mapped[i])))
            merged_diamondout = join(output_dir,
                    '{name}.singletons.nr.diamondout'.format(name=readset.name))

            jobs.append(Job(
                name='{step}.{name}.merge'.format(step=self.diamond_align_singletons.__name__,
                                                name=readset.name),
                input_files=[mapped[1], mapped[2]],
                output_files=[merged_diamondout],
                command='cat {out1} {out2} > {merged}'.format(out1=mapped[1],
                                                                out2=mapped[2],
                                                                merged=merged_diamondout)))


        return jobs

    def diamond_contigs_get_tophits(self):
        """
        Extract top his from contig DIAMOND serach output

        Input:
        {readset}.contigs.IDs_length.txt
        {readset}.contigs.nr.diamondout

        Output:
        {readset}.contigs.nr_diamond_pairs.txt
        {readset}.contigs.nr_diamond_IDs.txt
        """
        jobs = []

        input_prefix = 'contigs'
        output_prefix = 'contigs'

        for readset in self.readsets:
            input_dir = join(input_prefix, readset.name)
            output_dir = join(output_prefix, readset.name)

            length_file = join(input_dir,
                    '{name}.contigs.IDs_length.txt'.format(name=readset.name))
            diamondout = join(input_dir,
                    '{name}.contigs.nr.diamondout'.format(name=readset.name))

            tophits = join(output_dir,
                    '{name}.contigs.nr_diamond_pairs.txt'.format(name=readset.name))
            id_file = join(output_dir,
                    '{name}.contigs.nr_diamond_IDs.txt'.format(name=readset.name))

            jobs.append(Job(
                name='{step}.{name}.get_tophits'.format(step=self.diamond_contigs_get_tophits.__name__,
                                                        name=readset.name),
                input_files=[length_file, diamondout],
                output_files=[tophits, id_file],
                command='python {script_path}/get_blast_tophits.py '
                        '--length {length_file} '
                        '--infile {diamondout} '
                        '--outfile {tophits} '
                        '--id-file {id_file} '
                        '--cutoff-type 1 '
                        '--cutoff0 100 '
                        '--cutoff1 85 '
                        '--cutoff2 65 '
                        '--cutoff3 60 '
                        '--diamond true'.format(script_path=self.script_path,
                                                length_file=length_file,
                                                diamondout=diamondout,
                                                tophits=tophits,
                                                id_file=id_file)))

        return jobs

    def diamond_singletons_get_tophits(self):
        """
        Align singletons against non-redundant protein database using DIAMOND

        Input:
        {readset}.singletons_length.txt
        {readset}.singletons.nr.diamondout

        Output:
        {readset}.singletons.nr_sorted.diamondout
        {readset}.singletons.nr_diamond_pairs.txt
        {readset}.singletons.nr_diamond_IDs.txt
        """
        jobs = []

        input_prefix = 'contigs'
        output_prefix = 'contigs'

        for readset in self.readsets:
            input_dir = join(input_prefix, readset.name)
            output_dir = join(output_prefix, readset.name)

            length_file = join(input_dir,
                    '{name}.singletons_length.txt'.format(name=readset.name))
            diamondout = join(input_dir,
                    '{name}.singletons.nr.diamondout'.format(name=readset.name))
            sorted_out = join(output_dir,
                    '{name}.singletons.nr_sorted.diamondout'.format(name=readset.name))

            # Sort the diamondout file
            jobs.append(Job(
                name='{step}.{name}.sort'.format(step=self.diamond_singletons_get_tophits.__name__,
                                                    name=readset.name),
                input_files=[diamondout],
                output_files=[sorted_out],
                command='python {script_path}/sort_blastout.py '
                        '--infile {diamondout} '
                        '--sort {sorted_out} '
                        '--cutoff 10'.format(script_path=self.script_path,
                                                diamondout=diamondout,
                                                sorted_out=sorted_out)))
            # Extract tophits from the sorted file
            tophits = join(output_dir,
                    '{name}.singletons.nr_diamond_pairs.txt'.format(name=readset.name))
            id_file = join(output_dir,
                    '{name}.singletons.nr_diamond_IDs.txt'.format(name=readset.name))

            jobs.append(Job(
                name='{step}.{name}.get_tophits'.format(step=self.diamond_singletons_get_tophits.__name__,
                                                        name=readset.name),
                input_files=[length_file, sorted_out],
                output_files=[tophits, id_file],
                command='python {script_path}/get_blast_tophits.py '
                        '--length {length_file} '
                        '--infile {sorted_out} '
                        '--outfile {tophits} '
                        '--id-file {id_file} '
                        '--cutoff-type 1 '
                        '--cutoff0 100 '
                        '--cutoff1 85 '
                        '--cutoff2 65 '
                        '--cutoff3 60 '
                        '--diamond true'.format(script_path=self.script_path,
                                                length_file=length_file,
                                                sorted_out=sorted_out,
                                                tophits=tophits,
                                                id_file=id_file)))

        return jobs

    def generate_microbial_sequence(self):
        """
        Generate sequence file of mapped gene from BWA and BLAT mapping result

        Input:
        {readset}.contigs.micro_cds_id_bwa.txt
        {readset}.contigs.n_micro_cds_id_blat.txt
        {readset}.singletons.micro_cds_id_bwa.txt
        {readset}.singletons.n_micro_cds_id_blat.txt

        Output:
        {readset}.bwablat_merged_unique.txt
        {readset}.microbial_cds_sub.fasta
        {readset}.microbial_cds_sub_IDs_length.txt
        """
        jobs = [];

        input_prefix = 'contigs'
        output_prefix = 'contigs'

        for readset in self.readsets:
            input_dir = join(input_prefix, readset.name)
            output_dir = join(output_prefix, readset.name)

            id1 = join(input_dir,
                    '{name}.contigs.micro_cds_id_bwa.txt'.format(name=readset.name))
            id2 = join(input_dir,
                    '{name}.contigs.n_micro_cds_id_blat.txt'.format(name=readset.name))
            id3 = join(input_dir,
                    '{name}.singletons.micro_cds_id_bwa.txt'.format(name=readset.name))
            id4 = join(input_dir,
                    '{name}.singletons.n_micro_cds_id_blat.txt'.format(name=readset.name))

            unique_hit = join(output_dir,
                    '{name}.bwablat_merged_unique.txt'.format(name=readset.name))

            # Generate merged and unique hit id files
            jobs.append(Job(
                name='{step}.{name}.merge_id'.format(step=self.generate_microbial_sequence.__name__,
                                                        name=readset.name),
                input_files=[id1,id2,id3,id4],
                output_files=[unique_hit],
                command='python {script_path}/merge_id_files.py '
                        '--alignment bwa blat '
                        '--location {location} '
                        '--unique {unique_hit}'.format(script_path=self.script_path,
                                                            location=output_dir,
                                                            unique_hit=unique_hit)))

            # Generate sequence file using samtools
            blast_db = config.param(self.generate_microbial_sequence.__name__,
                    'blast_db', type='filepath')
            sequenced_fasta = join(output_dir,
                    '{name}.microbial_cds_sub.fasta'.format(name=readset.name))

            jobs.append(Job(
                name='{step}.{name}.generate_sequence'.format(step=self.generate_microbial_sequence.__name__,
                                                                name=readset.name),
                input_files=[unique_hit, blast_db],
                output_files=[sequenced_fasta],
                module_entries=[[self.generate_microbial_sequence.__name__,
                    'module_samtools']],
                command='xargs samtools faidx {blast} < {unique_hit} '
                        '> {output}'.format(blast=blast_db,
                                            unique_hit=unique_hit,
                                            output=sequenced_fasta)))

            # Generate length file
            length = join(output_dir,
                    '{name}.microbial_cds_sub_IDs_length.txt'.format(name=readset.name))

            jobs.append(Job(
                name='{step}.{name}.get_length'.format(step=self.generate_microbial_sequence.__name__,
                                                        name=readset.name),
                input_files=[sequenced_fasta],
                output_files=[length],
                command='python {script_path}/main_get_sequence_length.py '
                        '--fasta {input_fasta} '
                        '--output {length}'.format(script_path=self.script_path,
                                                    input_fasta=sequenced_fasta,
                                                    length=length)))

        return jobs

    def get_topbachit_contigs(self):
        """
        Extract top hits from multiple proteins with same score

        Input:
        {readset}.contigs.nr_diamond_IDs.txt
        {readset}.contigs.nr_diamond_pairs.txt

        Output:
        {readset}.contigs.nr_diamond_pairs_sub.txt
        {readset}.contigs.nr_diamond_hitsID_sub.txt
        {readset}.contigs.nr_diamond_hitsID_bacsub.txt
        """
        jobs = []

        input_prefix = 'contigs'
        output_prefix = 'contigs'

        for readset in self.readsets:
            input_dir = join(input_prefix, readset.name)
            output_dir = join(output_prefix, readset.name)

            blastdb = config.param(self.get_topbachit_contigs.__name__,
                    'blast_db', type='filepath')
            id_file = join(input_dir,
                    '{name}.contigs.nr_diamond_IDs.txt'.format(name=readset.name))
            pair_file = join(input_dir,
                    '{name}.contigs.nr_diamond_pairs.txt'.format(name=readset.name))

            pairs_sub = join(output_dir,
                    '{name}.contigs.nr_diamond_pairs_sub.txt'.format(name=readset.name))
            hitsID_sub = join(output_dir,
                    '{name}.contigs.nr_diamond_hitsID_sub.txt'.format(name=readset.name))
            hitsID_bacsub = join(output_dir,
                    '{name}.contigs.nr_diamond_hitsID_bacsub.txt'.format(name=readset.name))

            jobs.append(Job(
                name='{step}.{name}.1topbachits_contigs'.format(step=self.get_topbachit_contigs.__name__,
                                                                name=readset.name),
                input_files=[blastdb, id_file, pair_file],
                output_files=[pairs_sub, hitsID_bacsub, hitsID_sub],
                module_entries=[[self.get_topbachit_contigs.__name__,
                    'module_blast']],
                command='perl {script_path}/main_get_blast_fromfile_1topbachit.pl '
                        '{id_file} {pair_file} {hitsID_sub} {pairs_sub} '
                        '{hitsID_bacsub} {blastdb}'.format(script_path=self.script_path,
                                                            id_file=id_file,
                                                            pair_file=pair_file,
                                                            pairs_sub=pairs_sub,
                                                            hitsID_sub=hitsID_sub,
                                                            hitsID_bacsub=hitsID_bacsub,
                                                            blastdb=blastdb)))

        return jobs

    def get_topbachit_singletons(self):
        """
        Extract top hits from multiple proteins with same score

        Input:
        {readset}.singletons.nr_diamond_IDs.txt
        {readset}.singletons.nr_diamond_pairs.txt

        Output:
        {readset}.singletons.nr_diamond_pairs_sub.txt
        {readset}.singletons.nr_diamond_hitsID_sub.txt
        {readset}.singletons.nr_diamond_hitsID_bacsub.txt
        """
        jobs = []

        input_prefix = 'contigs'
        output_prefix = 'contigs'

        for readset in self.readsets:
            input_dir = join(input_prefix, readset.name)
            output_dir = join(output_prefix, readset.name)

            blastdb = config.param(self.get_topbachit_singletons.__name__,
                    'blast_db', type='filepath')
            id_file = join(input_dir,
                    '{name}.singletons.nr_diamond_IDs.txt'.format(name=readset.name))
            pair_file = join(input_dir,
                    '{name}.singletons.nr_diamond_pairs.txt'.format(name=readset.name))

            pairs_sub = join(output_dir,
                    '{name}.singletons.nr_diamond_pairs_sub.txt'.format(name=readset.name))
            hitsID_sub = join(output_dir,
                    '{name}.singletons.nr_diamond_hitsID_sub.txt'.format(name=readset.name))
            hitsID_bacsub = join(output_dir,
                    '{name}.singletons.nr_diamond_hitsID_bacsub.txt'.format(name=readset.name))

            jobs.append(Job(
                name='{step}.{name}.1topbachits_singletons'.format(step=self.get_topbachit_singletons.__name__,
                                                                name=readset.name),
                input_files=[blastdb, id_file, pair_file],
                output_files=[pairs_sub, hitsID_bacsub, hitsID_sub],
                module_entries=[[self.get_topbachit_singletons.__name__,
                    'module_blast']],
                command='perl {script_path}/main_get_blast_fromfile_1topbachit.pl '
                        '{id_file} {pair_file} {hitsID_sub} {pairs_sub} '
                        '{hitsID_bacsub} {blastdb}'.format(script_path=self.script_path,
                                                            id_file=id_file,
                                                            pair_file=pair_file,
                                                            pairs_sub=pairs_sub,
                                                            hitsID_sub=hitsID_sub,
                                                            hitsID_bacsub=hitsID_bacsub,
                                                            blastdb=blastdb)))

        return jobs

    def generate_nr_sequence(self):
        """
        Generate sequence file of mapped proteins from DIAMOND searches

        Input:
        {readset}.contigs.nr_diamond_hitsID_sub.txt
        {readset}.singletons.nr_diamond_hitsID_sub.txt
        {readset}.nr_diamond_hitsID_sub.txt

        Output:
        {readset}.nr_all_sub.fasta
        {readset}.nr_all_sub_IDs.txt
        {readset}.nr_all_sub_IDs_length.txt
        """
        jobs = []

        input_prefix = 'contigs'
        output_prefix = 'contigs'

        for readset in self.readsets:
            input_dir = join(input_prefix, readset.name)
            output_dir = join(output_prefix, readset.name)

            blastdb = config.param(self.generate_nr_sequence.__name__,
                    'blast_db', type='filepath')
            id_file1 = join(input_dir,
                    '{name}.contigs.nr_diamond_hitsID_sub.txt'.format(name=readset.name))
            id_file2 = join(input_dir,
                    '{name}.singletons.nr_diamond_hitsID_sub.txt'.format(name=readset.name))
            merged_id = join(input_dir,
                    '{name}.nr_diamond_hitsID_sub.txt'.format(name=readset.name))

            jobs.append(Job(
                name='{step}.{name}.merge_id'.format(step=self.generate_nr_sequence.__name__,
                                                        name=readset.name),
                input_files=[id_file1, id_file2],
                output_files=[merged_id],
                command='cat {id1} {id2} > {merged}'.format(id1=id_file1,
                                                            id2=id_file2,
                                                            merged=merged_id)))

            nr_fasta = join(output_dir,
                    '{name}.nr_all_sub.fasta'.format(name=readset.name))
            nr_id = join(output_dir,
                    '{name}.nr_all_sub_IDs.txt'.format(name=readset.name))

            jobs.append(Job(
                name='{step}.{name}.generate_nr_sub_fasta'.format(step=self.generate_nr_sequence.__name__,
                                                                    name=readset.name),
                input_files=[merged_id, blastdb],
                output_files=[nr_fasta, nr_id],
                module_entries=[[self.generate_nr_sequence.__name__,
                    'module_blast']],
                command='perl {script_path}/main_get_nr_sub.pl '
                        '{merged} {fasta} {out_id} '
                        '{blastdb}'.format(script_path=self.script_path,
                                            merged=merged_id,
                                            fasta=nr_fasta,
                                            out_id=nr_id,
                                            blastdb=blastdb)))

            # Get sequence length from generated fasta file
            length_file = join(output_dir,
                    '{name}.nr_all_sub_IDs_length.txt'.format(name=readset.name))

            jobs.append(Job(
                name='{step}.{name}.get_sequence_length'.format(step=self.generate_nr_sequence.__name__,
                                                                name=readset.name),
                input_files=[nr_fasta],
                output_files=[length_file],
                command='python {script_path}/main_get_sequence_length.py '
                        '--fasta {nr_fasta} '
                        '--output {length}'.format(script_path=self.script_path,
                                                    nr_fasta=nr_fasta,
                                                    length=length_file)))

        return jobs
#-------------------------------- Step 7 Ends Here ------------------------------------------


#-------------------------------- Step 8 Starts Here ------------------------------------------
    def align_genes_ecoli(self):
        """
        Match BWA and BLAT result to E. coli homologs

        Input:
        {readset}.microbial_cds_sub.fasta
        {readset}.microbial_cds_sub_IDs_length.txt

        Output:
        {readset}.microbial_cds_sub_ecoli_ppi.matches.daa
        {readset}.microbial_cds_sub_ecoli_ppi.diamondout
        {readset}.microbial_cds_sub_ecoli_ppi_pairs.txt
        {readset}.microbial_cds_sub_ecoli_ppi_IDs.txt
        """
        jobs = []

        input_prefix = 'contigs'
        output_prefix = 'contigs'

        for readset in self.readsets:
            input_dir = join(input_prefix, readset.name)
            output_dir = join(input_prefix, readset.name)
            tmp_dir = join(input_dir, 'diamond_tmp')
            tmp = join(tmp_dir, 'dir_created')

            ecoli_db = config.param(self.align_genes_ecoli.__name__,
                    'ecoli_db', type='filepath')

            input_fasta = join(input_dir,
                    '{name}.microbial_cds_sub.fasta'.format(name=readset.name))
            aligned = join(output_dir,
                    '{name}.microbial_cds_sub_ecoli_ppi.matches'.format(name=readset.name))
            aligned_renamed = join(output_dir,
                    '{name}.microbial_cds_sub_ecoli_ppi.matches.daa'.format(name=readset.name))

            # Align against E Coli database using Diamond
            jobs.append(Job(
                name='{step}.{name}.align'.format(step=self.align_genes_ecoli.__name__,
                                                    name=readset.name),
                input_files=[tmp, ecoli_db, input_fasta],
                output_files=[aligned_renamed],
                module_entries=[[self.align_genes_ecoli.__name__,
                    'module_diamond']],
                command='diamond blastx -p 8 -d {ecoli_db} -q {input_fasta} '
                        '-a {aligned} -t {tmp_dir} -e 10 -k 10'.format(
                                                                    ecoli_db=ecoli_db,
                                                                    input_fasta=input_fasta,
                                                                    aligned=aligned,
                                                                    tmp_dir=tmp_dir)))

            # Map read using diamond view
            reformatted = join(output_dir,
                    '{name}.microbial_cds_sub_ecoli_ppi.diamondout'.format(name=readset.name))

            jobs.append(Job(
                name='{step}.{name}.reformat_diamond'.format(step=self.align_genes_ecoli.__name__,
                                                                name=readset.name),
                input_files=[aligned_renamed],
                output_files=[reformatted],
                module_entries=[[self.align_genes_ecoli.__name__,
                    'module_diamond']],
                command='diamond view -a {aligned} -o {reformatted} '
                        '-f tab'.format(aligned=aligned_renamed,
                                        reformatted=reformatted)))

            # Extract top hits from diamond output
            length = join(input_dir,
                    '{name}.microbial_cds_sub_IDs_length.txt'.format(name=readset.name))
            tophit = join(output_dir,
                    '{name}.microbial_cds_sub_ecoli_ppi_pairs.txt'.format(name=readset.name))
            id_file = join(output_dir,
                    '{name}.microbial_cds_sub_ecoli_ppi_IDs.txt'.format(name=readset.name))

            jobs.append(Job(
                name='{step}.{name}.extract_tophits'.format(step=self.align_genes_ecoli.__name__,
                                                            name=readset.name),
                input_files=[reformatted],
                output_files=[tophit, id_file],
                command='python {script_path}/get_blast_1tophit.py '
                        '--length {length} '
                        '--infile {diamondout} '
                        '--outfile {tophit} '
                        '--id-file {id_file} '
                        '--cutoff-type 0 '
                        '--diamond true'.format(script_path=self.script_path,
                                                length=length,
                                                diamondout=reformatted,
                                                tophit=tophit,
                                                id_file=id_file)))

        return jobs

    def align_proteins_ecoli(self):
        """
        Match proteins identified through DIAMOND search to E. coli homologs

        Input:
        {readset}.nr_all_sub.fasta
        {readset}.nr_all_sub_IDs_length.txt
        Output:
        {readset}.nr_all_sub_ecoli_ppi.matches.daa
        {readset}.nr_all_sub_ecoli_ppi.diamondout
        {readset}.nr_all_sub_ecoli_ppi_pairs.txt
        {readset}.nr_all_sub_ecoli_ppi_IDs.txt
        """
        jobs = []

        input_prefix = 'contigs'
        output_prefix = 'contigs'

        for readset in self.readsets:
            input_dir = join(input_prefix, readset.name)
            output_dir = join(input_prefix, readset.name)
            tmp_dir = join(input_dir, 'diamond_tmp')
            tmp = join(tmp_dir, 'dir_created')

            ecoli_db = config.param(self.align_proteins_ecoli.__name__,
                    'ecoli_db', type='filepath')

            input_fasta = join(input_dir,
                    '{name}.nr_all_sub.fasta'.format(name=readset.name))
            aligned = join(output_dir,
                    '{name}.nr_all_sub_ecoli_ppi.matches'.format(name=readset.name))
            aligned_renamed = join(output_dir,
                    '{name}.nr_all_sub_ecoli_ppi.matches.daa'.format(name=readset.name))

            # Align against E Coli database using Diamond
            jobs.append(Job(
                name='{step}.{name}.align'.format(step=self.align_proteins_ecoli.__name__,
                                                    name=readset.name),
                input_files=[tmp, ecoli_db, input_fasta],
                output_files=[aligned_renamed],
                module_entries=[[self.align_proteins_ecoli.__name__,
                    'module_diamond']],
                command='diamond blastp -p 8 -d {ecoli_db} -q {input_fasta} '
                        '-a {aligned} -t {tmp_dir} -e 10 -k 10'.format(
                                                                    ecoli_db=ecoli_db,
                                                                    input_fasta=input_fasta,
                                                                    aligned=aligned,
                                                                    tmp_dir=tmp_dir)))

            # Map read using diamond view
            reformatted = join(output_dir,
                    '{name}.nr_all_sub_ecoli_ppi.diamondout'.format(name=readset.name))

            jobs.append(Job(
                name='{step}.{name}.reformat_diamond'.format(step=self.align_proteins_ecoli.__name__,
                                                                name=readset.name),
                input_files=[aligned_renamed],
                output_files=[reformatted],
                module_entries=[[self.align_proteins_ecoli.__name__,
                    'module_diamond']],
                command='diamond view -a {aligned} -o {reformatted} '
                        '-f tab'.format(aligned=aligned_renamed,
                                        reformatted=reformatted)))

            # Extract top hits from diamond output
            length = join(input_dir,
                    '{name}.nr_all_sub_IDs_length.txt'.format(name=readset.name))
            tophit = join(output_dir,
                    '{name}.nr_all_sub_ecoli_ppi_pairs.txt'.format(name=readset.name))
            id_file = join(output_dir,
                    '{name}.nr_all_sub_ecoli_ppi_IDs.txt'.format(name=readset.name))

            jobs.append(Job(
                name='{step}.{name}.extract_tophits'.format(step=self.align_proteins_ecoli.__name__,
                                                            name=readset.name),
                input_files=[reformatted],
                output_files=[tophit, id_file],
                command='python {script_path}/get_blast_1tophit.py '
                        '--length {length} '
                        '--infile {diamondout} '
                        '--outfile {tophit} '
                        '--id-file {id_file} '
                        '--cutoff-type 0 '
                        '--diamond true'.format(script_path=self.script_path,
                                                length=length,
                                                diamondout=reformatted,
                                                tophit=tophit,
                                                id_file=id_file)))

        return jobs

    def combine_ppi_results(self):
        """
        Combine mapped annotated genes and protein

        Input:
        {readset}.microbial_cds_sub_ecoli_ppi_pairs.txt
        {readset}.nr_all_sub_ecoli_ppi_pairs.txt
        Output:
        {readset}.PPI_pairs.txt
        """
        jobs = []

        input_prefix = 'contigs'
        output_prefix = 'contigs'

        for readset in self.readsets:
            input_dir = join(input_prefix, readset.name)
            output_dir = join(input_prefix, readset.name)

            gene_pair = join(input_dir,
                    '{name}.microbial_cds_sub_ecoli_ppi_pairs.txt'.format(name=readset.name))
            protein_pair = join(input_dir,
                    '{name}.nr_all_sub_ecoli_ppi_pairs.txt'.format(name=readset.name))
            combined = join(output_dir,
                    '{name}.PPI_pairs.txt'.format(name=readset.name))

            jobs.append(Job(
                name='{step}.{name}.combine_PPI_results'.format(step=self.combine_ppi_results.__name__,
                                                                name=readset.name),
                input_files=[gene_pair, protein_pair],
                output_files=[combined],
                command='perl {script_path}/main_combine_PPI_results.pl '
                        '{gene} {protein} {combined}'.format(script_path=self.script_path,
                                                            gene=gene_pair,
                                                            protein=protein_pair,
                                                            combined=combined)))


        return jobs
#-------------------------------- Step 8 Ends Here ------------------------------------------

#-------------------------------- Step 9 Starts Here ------------------------------------------
    def get_taxID_microbial(self):
        """
        Extract taxonomic information from generated microbial files

        Input
        {readset}.microbial_cds_sub_IDs_length.txt

        Output
        {readset}.microbial_cds_sub_IDs_map_taxid.txt

        """
        jobs = []

        input_prefix = 'contigs'
        output_prefix = 'contigs'

        for readset in self.readsets:
            input_dir = join(input_prefix, readset.name)
            output_dir = join(output_prefix, readset.name)

            # Map taxID to microbial_sub
            taxID = config.param(self.get_taxID_microbial.__name__, "taxID",
                    type='filepath')
            length = join(input_dir,
                    '{name}.microbial_cds_sub_IDs_length.txt'.format(name=readset.name))
            mapped = join(output_dir,
                    '{name}.microbial_cds_sub_IDs_map_taxid.txt'.format(name=readset.name))

            jobs.append(Job(
                name='{step}.{name}.map_taxID'.format(step=self.get_taxID_microbial.__name__,
                                                        name=readset.name),
                input_files=[taxID, length],
                output_files=[mapped],
                command='perl {script_path}/main_get_taxonID_microbial_cds.pl '
                        '{taxID} {length} {mapped}'.format(script_path=self.script_path,
                                                            taxID=taxID,
                                                            length=length,
                                                            mapped=mapped)))

        return jobs

    def get_taxID_nr(self):
        """
        Append taxID to E. coli mapped protein

        Input
        {readset}.nr_all_sub_IDs_length.txt
        {readset}.nr_all_sub_IDs.txt
        Output
        {readset}.nr_all_sub_IDs_map_taxid.txt
        {readset}.nr_all_sub_IDs_taxonID.txt
        """
        jobs = []

        input_prefix = 'contigs'
        output_prefix = 'contigs'

        for readset in self.readsets:
            input_dir = join(input_prefix, readset.name)
            output_dir = join(output_prefix, readset.name)

            blast_db = config.param(self.get_taxID_nr.__name__, 'blast_db',
                    type='filepath')
            length = join(input_dir,
                    '{name}.nr_all_sub_IDs_length.txt'.format(name=readset.name))
            id_file = join(input_dir,
                    '{name}.nr_all_sub_IDs.txt'.format(name=readset.name))
            mapped = join(output_dir,
                    '{name}.nr_all_sub_IDs_map_taxid.txt'.format(name=readset.name))
            taxonID = join(output_dir,
                    '{name}.nr_all_sub_IDs_taxonID.txt'.format(name=readset.name))

            jobs.append(Job(
                name='{step}.{name}.get_taxID'.format(step=self.get_taxID_nr.__name__,
                                                        name=readset.name),
                input_files=[length,id_file,blast_db],
                output_files=[mapped, taxonID],
                module_entries=[[self.get_taxID_nr.__name__, 'module_blast']],
                command='perl {script_path}/main_get_taxonID_nr.pl '
                        '{length} {id_file} {mapped} {taxonID} {blast_db}'.format(
                                                        script_path=self.script_path,
                                                        length=length,
                                                        id_file=id_file,
                                                        mapped=mapped,
                                                        taxonID=taxonID,
                                                        blast_db=blast_db)))

        return jobs

    def get_phylum_microbial(self):
        """
        Identify the phylum each annotated gene belongs to

        Input
        {readset}.microbial_cds_sub_IDs_map_taxid.txt

        Output
        {readset}.microbial_cds_sub_IDs_map_taxid_phylum.txt
        """
        jobs = []

        input_prefix = 'contigs'
        output_prefix = 'contigs'

        for readset in self.readsets:
            input_dir = join(input_prefix, readset.name)
            output_dir = join(output_prefix, readset.name)

            id_file = config.param(self.get_phylum_microbial.__name__,
                    'taxID_all', type='filepath')
            microbial_taxID = join(input_dir,
                    '{name}.microbial_cds_sub_IDs_map_taxid.txt'.format(name=readset.name))
            mapped = join(output_dir,
                    '{name}.microbial_cds_sub_IDs_map_taxid_phylum.txt'.format(name=readset.name))

            jobs.append(Job(
                name='{step}.{name}.get_phylum'.format(step=self.get_phylum_microbial.__name__,
                                                        name=readset.name),
                input_files=[id_file, microbial_taxID],
                output_files=[mapped],
                command='perl {script_path}/main_get_phylum.pl '
                        '{id_file} {microbial_taxID} {mapped}'.format(script_path=self.script_path,
                                                                        id_file=id_file,
                                                                        microbial_taxID=microbial_taxID,
                                                                        mapped=mapped)))

        return jobs

    def get_phylum_nr(self):
        """
        Identify the phylum each E. coli mapped protein belongs to

        Input
        {readset}.nr_all_sub_IDs_map_taxid.txt

        Output
        {readset}.nr_all_sub_IDs_map_taxid_phylum.txt
        """
        jobs = []

        input_prefix = 'contigs'
        output_prefix = 'contigs'

        for readset in self.readsets:
            input_dir = join(input_prefix, readset.name)
            output_dir = join(output_prefix, readset.name)

            id_file = config.param(self.get_phylum_nr.__name__,
                    'taxID_all', type='filepath')
            nr_taxID = join(input_dir,
                    '{name}.nr_all_sub_IDs_map_taxid.txt'.format(name=readset.name))
            mapped = join(output_dir,
                    '{name}.nr_all_sub_IDs_map_taxid_phylum.txt'.format(name=readset.name))

            jobs.append(Job(
                name='{step}.{name}.get_phylum'.format(step=self.get_phylum_nr.__name__,
                                                        name=readset.name),
                input_files=[id_file, nr_taxID],
                output_files=[mapped],
                command='perl {script_path}/main_get_phylum.pl '
                        '{id_file} {nr_taxID} {mapped}'.format(script_path=self.script_path,
                                                                id_file=id_file,
                                                                nr_taxID=nr_taxID,
                                                                mapped=mapped)))

        return jobs

    def get_mapped_geneIDs_microbial(self):
        """
        

        Input
        {readset}.microbial_cds_sub_IDs_length.txt'
        {readset}.contigs.IDs_length.txt
        {readset}.contigs.micro_cds_id_bwa.txt
        {readset}.singletons.micro_cds_id_bwa.txt
        {readset}.contigs.n_micro_cds_id_blat.txt
        {readset}.singletons.n_micro_cds_id_blat.txt

        Output
        {readset}.microbial_cds_sub_bwablat_pairs.txt
        {readset}.microbial_cds_sub_IDs_counts.txt

        """
        jobs = []

        input_prefix = 'contigs'
        output_prefix = 'contigs'

        for readset in self.readsets:
            input_dir = join(input_prefix, readset.name)
            output_dir = join(output_prefix, readset.name)

            microbialID = join(input_dir,
                    '{name}.microbial_cds_sub_IDs_length.txt'.format(name=readset.name))
            contigsID = join(input_dir,
                    '{name}.contigs.IDs_length.txt'.format(name=readset.name))

            contigs_bwa_pairs = join(input_dir,
                    '{name}.contigs.micro_cds_id_bwa.txt'.format(name=readset.name))
            singletons_bwa_pairs = join(input_dir,
                    '{name}.singletons.micro_cds_id_bwa.txt'.format(name=readset.name))
            contigs_blat_pairs = join(input_dir,
                    '{name}.contigs.n_micro_cds_id_blat.txt'.format(name=readset.name))
            singletons_blat_pairs = join(input_dir,
                    '{name}.singletons.n_micro_cds_id_blat.txt'.format(name=readset.name))

            combined_pair = join(output_dir,
                    '{name}.microbial_cds_sub_bwablat_pairs.txt'.format(name=readset.name))

            # Combine pair files
            jobs.append(Job(
                name='{step}.{name}.combine_pairs'.format(step=self.get_mapped_geneIDs_microbial.__name__,
                                                            name=readset.name),
                input_files=[contigs_bwa_pairs,singletons_bwa_pairs,contigs_blat_pairs,singletons_blat_pairs],
                output_files=[combined_pair],
                command="cat {contigs_bwa} {singletons_bwa} {contigs_blat} "
                        "{singletons_blat} | "
                        "sed -e '2,${{ /^@/d }}' > {combined}".format(
                                                contigs_bwa=contigs_bwa_pairs,
                                                singletons_bwa=singletons_bwa_pairs,
                                                contigs_blat=contigs_blat_pairs,
                                                singletons_blat=singletons_blat_pairs,
                                                combined=combined_pair)))

            # Map gene IDs
            mapped_gene = join(output_dir,
                    '{name}.microbial_cds_sub_IDs_counts.txt'.format(name=readset.name))

            jobs.append(Job(
                name='{step}.{name}.map_geneIDs'.format(step=self.get_mapped_geneIDs_microbial.__name__,
                                                        name=readset.name),
                input_files=[microbialID,contigsID,combined_pair],
                output_files=[mapped_gene],
                command='perl {script_path}/main_get_mapped_genesID_counts.pl '
                        '{microbialID} {contigsID} {combined_pair} '
                        '{mapped_gene}'.format(script_path=self.script_path,
                                                microbialID=microbialID,
                                                contigsID=contigsID,
                                                combined_pair=combined_pair,
                                                mapped_gene=mapped_gene)))

        return jobs

    def get_mapped_geneIDs_nr(self):
        """

        Input
        {readset}.nr_all_sub_IDs_length.txt
        {readset}.contigs.IDs_length.txt
        {readset}.contigs.nr_diamond_pairs_sub.txt
        {readset}.singletons.nr_diamond_pairs_sub.txt

        Output
        {readset}.nr_all_sub_combined_pairs_sub.txt
        {readset}.nr_all_sub_IDs_counts.txt

        """
        jobs = []

        input_prefix = 'contigs'
        output_prefix = 'contigs'

        for readset in self.readsets:
            input_dir = join(input_prefix, readset.name)
            output_dir = join(output_prefix, readset.name)

            nrID = join(input_dir,
                    '{name}.nr_all_sub_IDs_length.txt'.format(name=readset.name))
            contigsID = join(input_dir,
                    '{name}.contigs.IDs_length.txt'.format(name=readset.name))

            contigs_nr_pairs = join(input_dir,
                    '{name}.contigs.nr_diamond_pairs_sub.txt'.format(name=readset.name))
            singletons_nr_pairs = join(input_dir,
                    '{name}.singletons.nr_diamond_pairs_sub.txt'.format(name=readset.name))
            combined_pairs = join(output_dir,
                    '{name}.nr_all_sub_combined_pairs_sub.txt'.format(name=readset.name))

            # Combine pairs files
            jobs.append(Job(
                name='{step}.{name}.combine_pairs'.format(step=self.get_mapped_geneIDs_nr.__name__,
                                                            name=readset.name),
                input_files=[contigs_nr_pairs,singletons_nr_pairs],
                output_files=[combined_pairs],
                command='cat {contigs_nr} {singletons_nr} > {combined}'.format(
                                                                        contigs_nr=contigs_nr_pairs,
                                                                        singletons_nr=singletons_nr_pairs,
                                                                        combined=combined_pairs)))

            # Get mapped geneIDs
            mapped_gene = join(output_dir,
                    '{name}.nr_all_sub_IDs_counts.txt'.format(name=readset.name))

            jobs.append(Job(
                name='{step}.{name}.map_geneIDs'.format(step=self.get_mapped_geneIDs_nr.__name__,
                                                        name=readset.name),
                input_files=[nrID,contigsID,combined_pairs],
                output_files=[mapped_gene],
                command='perl {script_path}/main_get_mapped_genesID_counts.pl '
                        '{nrID} {contigsID} {combined_pairs} '
                        '{mapped_gene}'.format(script_path=self.script_path,
                                                nrID=nrID,
                                                contigsID=contigsID,
                                                combined_pairs=combined_pairs,
                                                mapped_gene=mapped_gene)))

        return jobs

    def get_mapped_gene_table_microbial(self):
        """

        Input
        {readset}.microbial_cds_sub_IDs_length.txt
        {readset}.microbial_cds_sub_IDs_map_taxid_phylum.txt
        {readset}.microbial_cds_sub_IDs_counts.txt
        {readset}.PPI_pairs.txt

        Output
        {readset}.microbial_cds_sub_table_counts.txt

        """
        jobs = []

        input_prefix = 'contigs'
        output_prefix = 'contigs'

        for readset in self.readsets:
            input_dir = join(input_prefix, readset.name)
            output_dir = join(output_prefix, readset.name)

            length = join(input_dir,
                    '{name}.microbial_cds_sub_IDs_length.txt'.format(name=readset.name))
            phylum = join(input_dir,
                    '{name}.microbial_cds_sub_IDs_map_taxid_phylum.txt'.format(name=readset.name))
            count = join(input_dir,
                    '{name}.microbial_cds_sub_IDs_counts.txt'.format(name=readset.name))
            ppi = join(input_dir,
                    '{name}.PPI_pairs.txt'.format(name=readset.name))
            table = join(output_dir,
                    '{name}.microbial_cds_sub_table_counts.txt'.format(name=readset.name))

            jobs.append(Job(
                name='{step}.{name}.get_mapped_gene_table_microbial'.format(step=self.get_mapped_gene_table_microbial.__name__,
                                                                            name=readset.name),
                input_files=[length,phylum,count,ppi],
                output_files=[table],
                command='perl {script_path}/main_get_mapped_gene_table.pl '
                        '{length} {phylum} {count} {ppi} {table}'.format(script_path=self.script_path,
                                                                        length=length,
                                                                        phylum=phylum,
                                                                        count=count,
                                                                        ppi=ppi,
                                                                        table=table)))

        return jobs

    def get_mapped_gene_table_nr(self):
        """

        Input
        {readset}.nr_all_sub_IDs_length.txt
        {readset}.nr_all_sub_IDs_map_taxid_phylum.txt
        {readset}.nr_all_sub_IDs_counts.txt
        {readset}.PPI_pairs.txt

        Output
        {readset}.nr_all_sub_table_counts.txt

        """
        jobs = []

        input_prefix = 'contigs'
        output_prefix = 'contigs'

        for readset in self.readsets:
            input_dir = join(input_prefix, readset.name)
            output_dir = join(output_prefix, readset.name)

            length = join(input_dir,
                    '{name}.nr_all_sub_IDs_length.txt'.format(name=readset.name))
            phylum = join(input_dir,
                    '{name}.nr_all_sub_IDs_map_taxid_phylum.txt'.format(name=readset.name))
            count = join(input_dir,
                    '{name}.nr_all_sub_IDs_counts.txt'.format(name=readset.name))
            ppi = join(input_dir,
                    '{name}.PPI_pairs.txt'.format(name=readset.name))
            table = join(output_dir,
                    '{name}.nr_all_sub_table_counts.txt'.format(name=readset.name))

            jobs.append(Job(
                name='{step}.{name}.get_mapped_gene_table_nr'.format(step=self.get_mapped_gene_table_nr.__name__,
                                                                        name=readset.name),
                input_files=[length,phylum,count,ppi],
                output_files=[table],
                command='perl {script_path}/main_get_mapped_gene_table.pl '
                        '{length} {phylum} {count} {ppi} {table}'.format(script_path=self.script_path,
                                                                        length=length,
                                                                        phylum=phylum,
                                                                        count=count,
                                                                        ppi=ppi,
                                                                        table=table)))

        return jobs

    def generate_RPKM(self):
        """
        Calculate a normalized expression value for each gene and protein

        Input
        {readset}.microbial_cds_sub_table_counts.txt
        {readset}.nr_all_sub_table_counts.txt

        Output
        {readset}.table_counts_all

        """
        jobs = []

        input_prefix = 'contigs'
        output_prefix = 'contigs'

        for readset in self.readsets:
            input_dir = join(input_prefix, readset.name)
            output_dir = join(output_prefix, readset.name)

            microbial_count = join(input_dir,
                    '{name}.microbial_cds_sub_table_counts.txt'.format(name=readset.name))
            nr_count = join(input_dir,
                    '{name}.nr_all_sub_table_counts.txt'.format(name=readset.name))
            combined_count = join(output_dir,
                    '{name}.table_counts_all'.format(name=readset.name))

            jobs.append(Job(
                name='{step}.{name}.combine_counts'.format(step=self.generate_RPKM.__name__,
                                                            name=readset.name),
                input_files=[microbial_count, nr_count],
                output_files=[combined_count],
                command='cat {microbial} {nr} > {combined}'.format(microbial=microbial_count,
                                                                    nr=nr_count,
                                                                    combined=combined_count)))

            RPKM = join(output_dir,
                    '{name}.table_RPKM_all.txt'.format(name=readset.name))

            jobs.append(Job(
                name='{step}.{name}.generate_RPKM'.format(step=self.generate_RPKM.__name__,
                                                            name=readset.name),
                input_files=[combined_count],
                output_files=[RPKM],
                command='perl {script_path}/main_get_mapped_gene_table_RPKM.pl '
                        '{combined} {RPKM}'.format(script_path=self.script_path,
                                                    combined=combined_count,
                                                    RPKM=RPKM)))

        return jobs
#-------------------------------- Step 9 Ends Here ------------------------------------------

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
            self.align_to_contigs,
            self.identify_contigs_reads,
            self.extract_singletons,  # 18
            self.get_mapping_table,
            self.bwa_align_contigs,
            self.bwa_identify_contigs,  # 21
            self.bwa_contigs_select_reads,
            self.bwa_align_singletons,
            self.bwa_identify_singletons,  # 24
            self.bwa_singletons_select_reads,
            self.blat_search_contigs,
            self.blat_search_singletons,  # 27
            self.process_contigs,
            self.process_singletons,
            self.diamond_align_contigs, # 30
            self.diamond_align_singletons,
            self.diamond_contigs_get_tophits,
            self.diamond_singletons_get_tophits,  # 33
            self.generate_microbial_sequence,
            self.get_topbachit_contigs,
            self.get_topbachit_singletons,  # 36
            self.generate_nr_sequence,
            self.align_genes_ecoli,
            self.align_proteins_ecoli,  # 39
            self.combine_ppi_results,
            self.get_taxID_microbial,
            self.get_taxID_nr,  # 42
            self.get_phylum_microbial,
            self.get_phylum_nr,
            self.get_mapped_geneIDs_microbial,  # 45
            self.get_mapped_geneIDs_nr,
            self.get_mapped_gene_table_microbial,
            self.get_mapped_gene_table_nr,  # 48
            self.generate_RPKM
        ]


if __name__ == '__main__':
    Metatranscriptomics()
