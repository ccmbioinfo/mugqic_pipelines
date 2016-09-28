#!/bin/bash

/hpf/projects/brudno/lmichael/rrbs_testrun/mugqic-2.2.0/pipelines/episeq/episeq.py \
 -s 1-6\
 -o /hpf/projects/brudno/lmichael/rrbs_testrun/output \
 -j pbs \
 -l debug \
 -d /hpf/projects/brudno/lmichael/rrbs_testrun/episeq.design \
 -r /hpf/projects/brudno/lmichael/rrbs_testrun/episeq.readset \
 -c episeq.ini 1> qsub.sh 2> debug.log 

chmod +x qsub.sh

