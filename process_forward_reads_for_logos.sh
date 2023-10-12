#!/bin/bash
module load cutadapt

input=$1
name=$(basename ${input%.fastq.gz})  
cutadapt  -n 2 -j 6 -m 26 \
--discard-untrimmed \
-g AGATGTGTATAAGAGACAG -a TGGATTGCGGGAAACGAG \
-o ./"${name}"_trimmed \
"${name}".fastq.gz
