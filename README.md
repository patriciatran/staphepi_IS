# *Staphylococcus epidermidis* project

# About this repository

This repository contains code to for the identification of insertion site elements in publicly available *Staph epidermidis (S. epi)* genomes.

The dataset was obtained from NCBI, downloaded on 06/04/2025, resulting in 225 genomes, using this exact command from the `ncbi-datasets` CLI tool. We only chose complete assemblies, no metagenome-assembled genomes (MAGs), and with the assembly source being RefSeq.

```         
datasets download genome taxon 'Staphylococcus epidermidis' --assembly-level complete --include genome,cds,protein,gbff --mag exclude --assembly-source 'RefSeq'
```

# Code overview

This repository contains for the bioinformatics analysis descripted in the Methods section "Identification of defense systems and insertion sites (IS) across S. epidermidis genomes" (Line 691-708), and discussed at lines 219-229 of the current version 1 of the preprint (see References below). The relevant figure in the Manuscript is Supplemental Figure 4.

# Repository folder and files description

This github repository contains:

- README.md : this file

Code files:

- Code.Rmd : Step-by-step description of the analyses involved, from downloading the raw data to running the python and R code to make the figures.

- extract_dnaA_to_1Mbp_downstream_rlmH.py : Python code to extract a relevant region of the genomes, using the raw data from NCBI.

Folders:

- blast/ : Query and output of the BLAST analysis for the IS131-meclike sequence. The database searched all.fasta is too large to upload to Github.

- figures/ : Figures generated from the R code, including extra figures not used in the paper.

- htcondor_scripts/ : HTcondor scripts to run ISEScan and DefenseFinder on CHTC using htcondor.

- NCBI-accessions.txt : List of the NCBI accession IDs downloaded using ncbi-datasets (raw input data)

- tables/: Contains a summarize version of the DefenseFinder output, and a table summarizing the coordinates of the dnaA and rlmH genes for the 1Mbp genome fragment.

# Analysis overview

## 1. Download input genomes and prepare the data

1.  Gather the genomes
2.  Reorganize all the genomes to have the dnaA gene, then whatever is between it and rlmH, then and extract the first 1Mbp (1,000,000bp) downstream of it.

## 2. Blast of IS431-meclike

3.  Use blastn to with the 431-meclike.fasta (taken from ISFinder website) against the reoriented genomes.
4.  Use R to plot the data.

## 3. ISEScan

Identify even more IS elements.

5.  Run `ISEScan` v.1.7.3 on the genomes to identify a variery of IS elements. Use HTCondor for high-throughput job submission.

6.  Use R to plot the data, only plot the ones that are complete (not partial)

## 4. DefenseFinder

7.  Use the same reorganized input genomes as were used for the ISEScan step above.

8.  Run `Defense-finder v.2.0.1` to identify the defense and anti-defense systems. Use HTCondor for high-throughput job submission.

9.  Use R to plot the data.

# Tools Referenced

- NCBI toolkit CLI: O’Leary et al., 2024 [publication](https://www.nature.com/articles/s41597-024-03571-y) ; [website](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/command-line-tools/download-and-install/)
- ISEScan: Xie et al., 2017 [publication](https://academic.oup.com/bioinformatics/article/33/21/3340/3930124)
- Defensefinder: Tesson et al., 2022 [publication](https://www.nature.com/articles/s41467-022-30269-9)

# Manuscript referenced

This is the code written by Patricia Tran for the project in collaboration with My Tran and Charlie Mo in the Department of Bacteriology at the University of Wisconsin-Madison.

Publication Pre-print:

My Tran, Angel Hernandez Viera, Jessica Rosu, Patricia Tran, Doran Goldman, Charlie Mo. *"Staphylococcus epidermidis* Genomic Plasticity Modulates Horizontal Gene Transfer" bioRxiv 2026.05.01.722213; doi: <https://doi.org/10.64898/2026.05.01.722213>
