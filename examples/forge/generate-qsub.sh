#!/bin/bash

/hpf/tools/centos6/mugqic-pipelines/latest_hpf/pipelines/forge/forge.py
 -s 1-33 \
 -o /hpf/tools/centos6/mugqic-pipelines/latest_hpf/examples/forge/output \
 -j pbs \
 -l debug \
 -r cheo.agilent.readset \
 -c cheo.ini 1> qsub.sh 2> debug.log

chmod +x qsub.sh

