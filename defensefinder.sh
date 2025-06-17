#!/bin/bash

SAMPLE="$1"

hmmsearch -h
defense-finder run -h

defense-finder run -o defense-finder_${SAMPLE} --models-dir /projects/bacteriology_tran_data/defense_finder_v2.0.0_models/ -a protein.faa

tar -cvf defense-finder_${SAMPLE}.tar.gz defense-finder_${SAMPLE}
