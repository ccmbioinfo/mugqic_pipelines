#!/bin/bash

# Total steps: 1-24
/hpf/tools/centos6/mugqic-pipelines/latest_hpf/pipelines/rnaseq/rnaseq.py \
 -s 1-23 \
 -o /hpf/tools/centos6/mugqic-pipelines/latest_hpf/examples/rnaseq/output \
 -j pbs \
 -l debug \
 -d rnaseq.arun.design \
 -r rnaseq.arun.readset \
 -c rnaseq.ini 1> qsub.sh 2> debug.log
# Add report: --report

chmod +x qsub.sh

