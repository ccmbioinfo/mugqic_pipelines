#!/usr/bin/env bash

module load seqtk

seqtk seq -a $1 > $2
