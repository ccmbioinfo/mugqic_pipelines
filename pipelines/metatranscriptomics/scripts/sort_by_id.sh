#!/usr/bin/env bash
# Sort fastq file $1 by its id
# Output fastq file on stdout

tmp=`mktemp`
cat $1 | paste - - - - | sort -k1,1 -S 3G | tr '\t' '\n' > $tmp
mv $tmp $1
