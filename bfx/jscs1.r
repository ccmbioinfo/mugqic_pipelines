#!/usr/bin/env R

library(JunctionSeq);
library(JctSeqData)

decoder <- read.table("jctseq/jctseq.design", header=TRUE, stringsAsFactors=FALSE);

gff.file <- "jctseq/jctseq.gff.gz"
