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

#### Understanding the table:

---
#### Installation:



#### Test data
