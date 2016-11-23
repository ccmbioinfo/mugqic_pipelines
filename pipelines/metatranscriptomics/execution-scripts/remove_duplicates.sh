#!/bin/bash

module load mugqic-pipelines/2.2.0

mkdir remove_duplicates

python /hpf/largeprojects/ccmbio/nreinhardt/mugqic_pipelines/pipelines/metatranscriptomics/scripts/main_get_derepli_IDs.py
