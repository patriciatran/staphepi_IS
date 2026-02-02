#!/bin/bash

SAMPLE="$1"
CPU="$2"

isescan.py --seqfile ${SAMPLE}.fasta --output results_${SAMPLE} --nthread ${CPU}

tar cvf results_${SAMPLE}.tar.gz results_${SAMPLE}

ls -lht
