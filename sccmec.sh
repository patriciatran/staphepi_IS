#!/bin/bash

SAMPLE="$1"

wget https://github.com/rpetit3/sccmec/raw/refs/heads/main/data/sccmec-targets.fasta
wget https://github.com/rpetit3/sccmec/raw/refs/heads/main/data/sccmec-regions.fasta
wget https://github.com/rpetit3/sccmec/raw/refs/heads/main/data/sccmec-targets.yaml
wget https://github.com/rpetit3/sccmec/raw/refs/heads/main/data/sccmec-regions.yaml


sccmec-main --input ${SAMPLE}.fna \
     -yt sccmec-targets.yaml \
     -yr sccmec-regions.yaml \
     -t sccmec-targets.fasta \
     -r sccmec-regions.fasta

tar -cf ${SAMPLE}_sccmec_output.tar.gz *.tsv