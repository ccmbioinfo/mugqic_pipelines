#!/usr/bin/env python

# Python Standard Modules
import os
import ntpath

# MUGQIC Modules
from core.config import *
from core.job import *

def compress_bam(input_bam, output_cram):

  return Job(
      [input_bam],
      [output_cram],
      [
          ['cramtools_compress_bam','module_java'],
          ['cramtools_compress_bam','module_cramtools']
      ],
      command="""\
cp {input_bam} $TMPDIR && \\
cd $TMPDIR  && \\
java -Xmx{ram} -jar $CRAMTOOLS_JAR cram \\
  --input-bam-file {tmp_input_bam} \\
  --output-cram-file {tmp_output_cram} \\
  --reference-fasta-file {reference_fasta} \\
  --capture-all-tags \\
  --lossy-quality-score-spec \*8 && \\
cp {tmp_output_cram} $PBS_O_INITDIR/{output_dir}/ && \\
cd $PBS_O_INITDIR""".format(
			ram = config.param('cramtools_compress_bam','ram'),
			input_bam = input_bam,
			output_cram = output_cram,
			reference_fasta = config.param('cramtools_compress_bam','genome_fasta', required = True),
      tmp_input_bam = "$TMPDIR/" + ntpath.split(input_bam)[1],
      tmp_output_cram = "$TMPDIR/" + ntpath.split(output_cram)[1],
      output_dir = ntpath.split(output_cram)[0]
		)
	)
