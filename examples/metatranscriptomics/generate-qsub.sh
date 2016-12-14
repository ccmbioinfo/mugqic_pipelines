#!/bin/bash

module load mugqic-pipelines/2.2.0
export PYTHONPATH=/hpf/largeprojects/ccmbio/nreinhardt/metatranscriptomics-test:$PYTHONPATH

script=../../pipelines/metatranscriptomics/metatranscriptomics.py 
readset=cow.readset 
config=config.ini 

output=output

python $script \
 -o $output \
 -r $readset\
 -s 1-15 \
 -j pbs \
 -l debug \
 -c $config \
 1> qsub.sh \
 --force 

chmod +x qsub.sh

