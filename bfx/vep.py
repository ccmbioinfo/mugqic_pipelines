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

#!/usr/bin/env python

# Python Standard Modules

# MUGQIC Modules
from core.config import config
from core.job import Job


def annotate(vcf, out_file):
    return Job(
        [vcf],
        [out_file],
        [
            ['DEFAULT', 'module_perl']
        ],
        command="""\
perl /hpf/tools/centos6/vep/82/scripts/variant_effect_predictor/variant_effect_predictor.pl --offline --dir_cache /hpf/tools/centos6/vep/source/cache/ --assembly GRCh37 --vcf --sift b --polyphen b --symbol --numbers --biotype --total_length --canonical --ccds --fields Consequence,Codons,Amino_acids,Gene,SYMBOL,Feature,EXON,PolyPhen,SIFT,Protein_position,BIOTYPE,CANONICAL,CCDS,RadialSVM_score,RadialSVM_pred,LR_score,LR_pred,CADD_raw,CADD_phred,Reliability_index,LoF,LoF_filter,LoF_flags --fork {processors} --force_overwrite --out_file STDOUT --stats_file {out_file}.summary.html| grep -v -- "- INFO: Disabling"
""".format(processors=16, out_file=out_file)
    )
    # TODO make configurable
