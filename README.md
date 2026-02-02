# Staph epi project


## About this repository

This repository contains code to for the identification of insertion site elements in publically available Staph epi genomes.

The dataset was obtained from NCBI, downloaded on 06/04/2025, resulting in 225 genomes, using this exact command from the `ncbi-datasets` CLI tool. 
We only chose complete assemblies, no metagenome-assembled genomes (MAGs), and with the assembly source being RefSeq. 

```
datasets download genome taxon 'Staphylococcus epidermidis' --assembly-level complete --include genome,cds,protein,gbff --mag exclude --assembly-source 'RefSeq'
```

## Ref
This is the code written by Patricia Tran for the project in. collaboration with My Tran and Charlie Mo in the Department of Bacteriology at the University of Wisconsin-Madison.

Publication: [to be announced]
