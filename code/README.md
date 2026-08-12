# IS elements:

## Prepare the data

1. Gather the genomes
2. Reorganize all the genomes to have the dnaA gene, then whatever is between it and rlmH, then and extract the first 1Mbp (1,000,000bp) downstream of it.

### Blast of IS431-meclike
3. Use blastn to with the 431-meclike.fasta (taken from ISFinder website) against the reoriented genome. 
4. Use R to plot the data.

### ISEScan
Identify even more IS elements.

5.  Run `ISEScan` v.1.7.3 on the genomes to identify a variery of IS elements. Use HTCondor for high-throughput job submission.

6. Use R to plot the data, only plot the ones that are complete (not partial)

# DefenseFinder

7. Use the same reorganized input genomes as were used for the ISEScan step above.

8. Run `Defense-finder v.2.0.1` to identify the defense and anti-defense systems. Use HTCondor for high-throughput job submission.

9. Use R to plot the data. 

### Reference:

- NCBI toolkit CLI: O’Leary et al., 2024 [publication](https://www.nature.com/articles/s41597-024-03571-y) ; [website](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/command-line-tools/download-and-install/)
- ISEScan: Xie et al., 2017 [publication](https://academic.oup.com/bioinformatics/article/33/21/3340/3930124)
- Defensefinder: Tesson et al., 2022 [publication](https://www.nature.com/articles/s41467-022-30269-9)
