#!/usr/bin/env python

# Python Standard Modules

# MUGQIC Modules
from core.config import *
from core.job import *

def compute_tdf(input, output):
    return Job(
        [input],
        [output],
        [
            ['compute_tdf', 'module_java'],
            ['compute_tdf', 'module_igvtools']
        ],
        command="""\
igvtools count -f min,max,mean \\
  {input} \\
  {output} \\
  {genome}""".format(
        input=input,
        output=output,
        genome=config.param('compute_tdf', 'igv_genome', type='filepath')
        )
    )
