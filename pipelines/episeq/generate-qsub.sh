#!/bin/bash

# Setup environment to run MUGQIC pipelines
module load mugqic-pipelines/2.2.0

# Setup path variable for python. Modify these paths to point the
export PYTHONPATH=/hpf/projects/brudno/lmichael/test2/:$PYTHONPATH

./episeq/episeq.py \
 -s 1-13 \
 -o output \
 -j pbs \
 -l debug \
 -d ./episeq.design \
 -r ./episeq.readset \
 -c episeq.ini 1> qsub.sh 2> debug.log

# Other flags
#--report || --clean

chmod +x qsub.sh

