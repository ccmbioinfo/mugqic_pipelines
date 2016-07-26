#!/bin/bash

# Total Steps: 1-29 \
/hpf/tools/centos6/mugqic-pipelines/latest_hpf/pipelines/dnaseq/dnaseq.py \
 -s 1-28 \
 -o /hpf/tools/centos6/mugqic-pipelines/latest_hpf/examples/dnaseq/output \
 -j pbs \
 -l debug \
 -r dnaseq.readset \
 -c dnaseq.hpf.ini 1> qsub.sh 2> debug.log

chmod +x qsub.sh

