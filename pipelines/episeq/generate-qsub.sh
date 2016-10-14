#!/bin/bash

/hpf/projects/brudno/lmichael/mugqic_pipelines/pipelines/episeq/episeq.py \
 -s 1-7 \
 -o /hpf/projects/brudno/lmichael/episeq_run/output \
 -j pbs \
 -l debug \
 -d /hpf/projects/brudno/lmichael/episeq_run/episeq.design \
 -r /hpf/projects/brudno/lmichael/episeq_run/episeq.readset \
 -c episeq.ini 1> qsub.sh 2> debug.log 

chmod +x qsub.sh

