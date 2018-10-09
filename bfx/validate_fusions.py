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


#validated_fusions = open("integrate-validated-fusions.txt", 'r')
#output_fusions = open("merged.cff.reann.dnasupp.bwafilter.30.cluster", 'r')


# determines which of the validate fusions the pipeline is able to detect using corresponding sample sequence data
# containing the corresponding fusions
def validate_fusions(in_dir, out_dir, validated_fusions_file, ini_section='fusion_stats', repeat_filter_section='repeat_filter'):

    # load seq_len used in repeat_filter step
    seq_len = config.param(repeat_filter_section, 'seq_len', type='int')

    bwafilter_cluster_file = os.path.join(in_dir, "merged.cff.reann.dnasupp.bwafilter." + str(seq_len) + ".cluster" )

    output_file = os.path.join(out_dir, "validate_fusion_stats.txt")

    return Job(
        [bwafilter_cluster_file],
        [output_file],
        [["merge_and_reannotate_cff_fusion", "module_fusiontools"]],
        command="""\
validate_fusion_stats.py \\
  {bwafilter_cluster_file} \\
  {validated_fusions_file} \\
  > {output_file}""".format(
        bwafilter_cluster_file=bwafilter_cluster_file,
        validated_fusions_file=validated_fusions_file,
        output_file=output_file,
        ),
        removable_files=[]
    )
