#!/usr/bin/env python

################################################################################
# Copyright (C) 2014, 2015 GenAP, McGill University and Genome Quebec Innovation Centre
#
# This file is part of MUGQIC Pipelines.
#
# MUGQIC Pipelines is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# MUGQIC Pipelines is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with MUGQIC Pipelines.  If not, see <http://www.gnu.org/licenses/>.
################################################################################

# Python Standard Modules
import os

# MUGQIC Modules
from core.config import *
from core.job import *

# generates statistics on the category/cluster file, generated my the merge_and_reannotate_cff_fusion step 

def fusion_stats(in_dir, out_dir, sampleinfo_file,  ini_section='fusion_stats', repeat_filter_section='repeat_filter'):

    # TODO: implement other_options flags in category_fusion_stats.py
    other_options = config.param(ini_section, 'other_options', required=False)
    
    # load seq_len used in repeat_filter step
    seq_len = config.param(repeat_filter_section, 'seq_len', type='int')

    #cluster_file = os.path.join(out_dir, "merged.cff.reann.dnasupp.cluster")
    bwafilter_cluster_file = os.path.join(in_dir, "merged.cff.reann.dnasupp.bwafilter." + str(seq_len) + ".cluster" )
    
    output_file = os.path.join(out_dir, "fusion_stats.txt") 

    return Job(
        [bwafilter_cluster_file],
        [output_file],
        [["merge_and_reannotate_cff_fusion", "module_fusiontools"]],
        command="""\
category_fusion_stats.py \\
  {bwafilter_cluster_file} \\
  {sampleinfo_file} \\
  > {out_file}""".format(
        bwafilter_cluster_file=bwafilter_cluster_file,
        sampleinfo_file=sampleinfo_file,
        out_file=output_file,
        ),
        removable_files=[]
    )

# category_fusion_stats.py testfiles/merged.cff.reann.dnasupp.bwafilter.30.cluster testfiles/sampleinfo

# cluster type, head gene, tail gene, max split read cnt, max spanning read cnt,  sample type, disease
#tools, inferred fusion category, head breakpoint on boundary, head breakpoint close to boundary, tail breakpoint on boundary,
#tail breakpoint close to boundary, dna support, samples

#1.  cluster type
#2.  head gene
#3.  tail gene
#4.  max split read cnt
#5.  max spanning read cnt
#6.  sample type
#7.  disease
#8.  tools
#9.  inferred fusion category
#10. head breakpoint on boundary
#11. head breakpoint close to boundary
#12. tail breakpoint on boundary
#13. tail breakpoint close to boundary
#14. dna support
#15. samples

