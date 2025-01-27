This is the simulation procedure to 

1. evaluate the accuracy of GND and AAD estimates ("codon" simulations).
2. evaluate the accuracy of the LGT-free model ("nt" simulations)

### Goal: 
* Generate a set of mock genomes by introducing genomic changes only through point mutations (no genome rearrangements or HGT).  
* Given a genome X, add mutations to it in either the AA space and back-translated them to nucleotides ("codon") or in nucleotide space directly ("nt") to get a mock genome X'. Vary the parameters of the evolutionary model to get multiple mock genomes.


### Evolutionary model: 
We use two different evolutionary models to generate the data, a codon-based model and a nucleotide-based model. The inputs to our model are:

 1. `p`: refers to AAD for codon based model and GND for nuleotide model,
 2. `alpha`: the alpha parameter for the gamma distribution,
 3. `model`: "codon" for codon-based model and "nt" for nucleotide-based mode.
 4. genome id, and 
 5. base directory, for the codon-based model.

To generate a mock genome `X'` that has `AAD(X,X') = p` using the codon-based model, we follow the following steps:

  1. use a gamma distribution (alpha = beta) to draw the relative rate for each gene. 
  
  2. sample without replacement nmus = p\*L amino acids in X, where L is the length of X. Each amino acid in X is selected with a probability determined by the rate of the gene it belongs to. Mutate these amino acids using the BLOSUM62 model to obtain X'.  
  
  3. back translate X' from amino acids to nucleotides. 

The generation process of the nucleotide-based model is similar to the codon-based model, with the difference that the mutations are introduced at the nucleotide level and then translated to codons.

We vary `p` to create a set of `X'` with different distances to `X`.

### Simulated data
1. Codon Simulations:
We selected 10 different base genomes (5 pairs), and for each base genome simulated 35 mock genomes from the reference (base genome) in the codon mode (`model = "codon"`) with different AADs from `AAD = 0.01` to `AAD=0.69` with step size `0.02`. For these simulations we used `alpha = 22` for the gamma distribution. The simulated genomes are stored in ```genome_pairs``` folder.
2.  Nucleotide Simulations:
We selected 29 different base genomes, and for each base genome we simulated a set of 31 mock genomes from the reference genome in the nucleotide model (`model = "nt"`) with different GNDs ranging from `GND=0` to `GND=0.20` with `alpha=5.28` for the gamma distribution. The simulated genomes will be stored in ```simulated_genomes``` folder.

 
Inside each folder, there are 3 files: 
  
  1. ```gnd_aad.txt``` - the GND and AAD distance of this mock genome to the reference genome.
  2. ```mutated.fasta``` - the simulated mock genome in DNA, and 
  3. ```mutated_genes.faa``` - the simulated mock genome in amino acid.

### Simulation procedure
To reproduce the "codon" simulations data, run ```evolve_codon_multi_level.sh```. This bash script will call ```evolve.py``` (with `model = "codon"` and `alpha = 22`) multiple times to sequentially simulate 35 mock genomes for each base genome.

To reproduce the "nucleotide" simulations data, run ```evolve_nt_multi_level.sh```. This bash script will call ```evolve.py``` (with `model = "nt"` and `alpha = 5.28`) multiple times to sequentially simulate 31 mock genomes for each base genome.

### Side notes  
  * Back translation: when there are multiple codons for an amino acid, if there is a codon that is identical to the original genome, we choose it; otherwise, we randomly select a codon with weight proportional to 4-d where d is its Hamming distance to the original genome. 
  * Non-mutated regions: we don't mutate (1) the start/end codons, (2) intergenic regions (~5% of the genome), and (3) overlapping regions of the genes - when encounter a pair of genes that overlap each other (possibly with different reading frames), we don't mutate the overlapping region
