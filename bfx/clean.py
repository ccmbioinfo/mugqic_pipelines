#!/usr/bin/env python

# Python Standard Modules
import os
import ntpath

# MUGQIC Modules
from core.config import *
from core.job import *

def remove_intermediate_files(input_crams):

  return Job(
      input_crams,
      [],
      [],
      command="""\
mkdir -p norRNA_crams && \\
mv {norRNA_crams} norRNA_crams/ &&
rm -rf alignment_1stPass && \\
rm -rf metrics/rnaseqRep && \\
rm -rf reference.Merged && \\
rm -rf trim""".format(
			norRNA_crams = " \\\n".join(str(cram) for cram in input_crams)
		)
	)

#rm -rf alignment && \\
