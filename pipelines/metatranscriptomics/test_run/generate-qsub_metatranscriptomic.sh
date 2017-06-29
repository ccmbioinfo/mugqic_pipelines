#!/bin/bash

/hpf/projects/brudno/daniel/mugqic_pipelines/pipelines/metatranscriptomics/metatranscriptomics.py \
 -s 1-49 \
 -o /hpf/projects/brudno/daniel/mugqic_pipelines/pipelines/metatranscriptomics/output \
 -j pbs \
 -l debug \
 -r cow.readset \
 -c config.ini 1> qsub_metatranscriptomic.sh #2>debug.log

chmod +x qsub_metatranscriptomic.sh

