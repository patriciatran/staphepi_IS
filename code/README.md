# IS elements:

## Prepare the data

1. Gather the genomes
2. Reorganize all the genomes to have the dnaA gene as the first gene on the genome, and extract the first 300,000bp downstream of it.

### Blast of IS431-meclike
3. Use blastn to with the 431-meclike.fasta against the reoriented genome. 
4. Use R to plot the data.

### ISEScan
Identify even more IS elements.

5.  Run ISEScan (https://github.com/xiezhq/ISEScan) on the genomes to identify a variery of IS elements. Use HTCondor for high-throughput job submission.

6. Use R to plot the data, only plot the ones that are complete (not partial)

# DefenseFinder

7. Use the same reorganized input genomes as were used for the ISEScan step above.

8. Run `Defense-finder v.2.0.1` to identify the defense and anti-defense systems. Use HTCondor for high-throughput job submission.

9. Use R to plot the data. 

