## Comparing shared ORFs to ANI per genome (on LGT-free simulated genomes)

---  
#### What this pipeline does: Tabulates the number of orthologous genes per genome set at different pct_id thresholds.
#### Why? To track how homology (or its detectability) erodes as genomes diverge in sequence space.

### Context within the study:
* This pipeline takes simulated LGT-free genomes made in the [GNDSim](GNDSim) section.
* It tabulates differences between the genomes (i.e. shared ORFs vs ANI).
* These tables serve as inputs for plots generated in the [GNDModel](GNDModel) section.
* In the study, we iterated this pipeline over 56 sets of simulated LGT-free genomes.

#### How it works:
1. Takes a set of simulated LGT-free genomes. (All simulations are copies of the same real genome but with different numbers of point-mutations added in silico).
2. Runs pyANI comparison of all genomes in the set to calculate Average Nucleotide Identity.
4. Uses prokka to call genes (i.e. ORFs) for each genome.
2. Creates a BLAST database of all ORFs in the genome set.
3. Runs a BLASTn query of EACH individual genome's ORFs against ALL genomes' ORFs in the set.
4. Summarizes each BLAST file per query genome, writes each summary per SAG to a file.
5. Combines all summaries into one table, combines with pyANI output. Full outputs [here](https://github.com/smirarab/GORG-LGT/tree/master/GNDModel/simulations)

#### Inputs and outputs:
* Inputs
  * This pipeline takes (at least 1) zipped set of simulated genomes. We analyzed 56 sets, found here [input](https://github.com/smirarab/GORG-LGT/tree/master/ORFvANI/input)
  * Expects each set as a tarball named something like "AG-359-G18_a22.tar.gz"
  * The tarball contains genomes named like "AG-899-G06_a22_gnd001.fasta", "AG-899-G06_a22_gnd002.fasta", etc.
  * What do the names mean?
    * The rootname [AG-359-G18] is the name of a real bacterial genome from the GORG-tropics dataset.
    * The first suffix (either [a5] or [a22]) represents the alpha value used in the model that simulated point mutations.
    * The last suffix (e.g. [gnd001] vs. [gnd002]) represents the simulated genome's GND from the real genome.
* Notes
  * Expects genome sequence names in default SPAdes format, like ">AG-359-G18_NODE_1"
* Outputs
  * For each input genome set, it produces a separate table. We generated 56 tables, found [here](https://github.com/smirarab/GORG-LGT/tree/master/GNDModel/simulations)

#### Understanding the output table:
What do the columns mean?

* <ins>qsag</ins>: the query genome used in BLASTn, e.g. "AG-359-G18_a22_gnd00000"
* <ins>ssag</ins>: the subject genome, e.g. "AG-359-G18_a22_gnd00005"
* <ins>qgene_count</ins>: number of ORFs (predicted by prokka) in query genome
* <ins>sgene_count</ins>: number of ORFs (predicted by prokka) in subject genome
* <ins>qgene_count_500-1500bp</ins>: number of ORFs for size range 500-1500 bp in query genome
* <ins>sgene_count_500-1500bp</ins>: number of ORFs for size range 500-1500 bp in subject genome
* <ins>total_hits</ins>: number of ORFs hitting between the query and subject genome (via BLASTn)
* <ins>mean_pident</ins>: mean % identity across all ORFs hitting from query to subject genome
* <ins>median_pident</ins>: median % identity of all " " "
* <ins>stdev_pident</ins>: standard deviation of all " " " "
* <ins>99.9_pid_orthologs</ins>: number of BLASTn hits with pid => 99.9 between query and subject genomes
* <ins>99.9_pid_500-1500bp_orthologs</ins>: number of " " " " for size range 500-1500
  * Note: See also columns for percent identities of 99.8, 99.7, 99.6, 99.5, 97.5 and 97
* <ins>ani</ins>: Average Nucleotide Identity between query and subject genome (via pyANI)
* <ins>ani_aln_cov_ab</ins>: Percent of query genome that pyANI could align to subject genome
* <ins>ani_aln_cov_ba</ins>: Percent of subject genome that pyANI could align to query genome
* <ins>GND</ins>: Genome-wide Nucleotide Difference between the subject and query genome (i.e. the inverse of their Average Nucleotide Identity)

---
#### Installation:
Install [Nextflow](https://www.nextflow.io/docs/latest/install.html) along with the a container manager (either Singularity or Docker).

Make a [nextflow.config](https://www.nextflow.io/docs/latest/config.html) file with configurations appropriate for your system.

#### To run this pipeline:
Create a working directory with the following items:
1. The [orfs_vs_ani.nf](https://github.com/smirarab/GORG-LGT/blob/master/ORFvANI/orfs_vs_ani.nf) nextflow script
2. The directory [input/](https://github.com/smirarab/GORG-LGT/tree/master/ORFvANI/input) containing the zipped, simulated genomes to compare
3. Your nextflow.config file

Then, run the pipeline

``nextflow run orfs_vs_ani.nf``


Upon running, nextflow will use the container manager (either Singularity or Docker--whichever you installed) to download and install the following dependencies:

pyANI version [0.2.9](https://quay.io/repository/biocontainers/pyani?tab=tags)

Prokka version [1.14.6](https://quay.io/repository/biocontainers/prokka?tab=tags)

BLAST version [2.7.1](https://quay.io/repository/biocontainers/blast)

