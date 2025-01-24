* [aad_distances.txt.xz](aad_distances.txt.xz): csv file (xz compressed) showing pairwise AAD between genomes

* [aad_tree.nwk](aad_tree.nwk): Results of hierarchical single linkage cluster performed on AAD using [TreeN93](https://github.com/niemasd/TreeN93/tree/a4e2bfc8a0bd573d484165b7e99a53aa5eb443b9). The missing AAD pairs are set to 60 (i.e. a very high AAD distance).
  ~~~bash
  xz -d aad_distances.txt.xz
  TreeN93.py -i aad_distances.txt -o aad_tree.nwk -m 60
  ~~~

* [astral.tre](astral.tre): ASTRAL-Pro tree created from gene trees available [here](../genetreeAAD/all-gene-trees.zip). 

* [RAxML.T16S-SINA-1000-bootcutoff_MVrooted](RAxML.T16S-SINA-1000-bootcutoff_MVrooted): The 16S gene tree computed using RAxML with 1000 bootstrapping replicates

* [RAxML_bootstrap.T16S-SINA-1000-bootcutoff](RAxML_bootstrap.T16S-SINA-1000-bootcutoff): 1000 bootstrapping replicates of the 16S gene tree computed using RAxML

* `RAxML.T16S-SINA-1000-bootcutoff_MVrooted_pruned`, 
  `RAxML.T23S-SINA-1000-bootcutoff_MVrooted_pruned`, 
  `aad_tree_pruned.nwk`, 
  `astral_pruned.tre`: pruned version of trees above; restrict the leafset of all trees (astral, S16, S23, aad_tree) to their intersection

* [aad_binning.py](aad_binning.py): the sampling strategy using the hierarchical single linkage clustering. Creates these files:

  * [16s_sampledtrees.zip](16s_sampledtrees.zip): 16S tree subsampled to various ranges (e.g., `16SS_037_052.nwk`), ASTRAL restricted to the same leaves (`astral_037_052.nwk`), and mean pdistance (i.e. AAD distances) of that subset (`pdist_037_052.txt`).

* [combined_data.txt](combined_data.txt): A file that contains:
  * mean pairwise AAD (`pdist`) for a subtree
  * quartet distance between the species tree and the gene tree (`qdistAstral`, given for both 16S and 23S)
    * Computed using [tqDist](https://users-cs.au.dk/cstorm/software/tqdist/) as adopted by [tripVote]()
 
    ~~~bash
    # to get the ASTRAL quartet score
     quartet_dist -v [pruned astal tree] [pruned 16S tree] | awk '{print $4;}
    ~~~
  * mean quartet distance btween bootstrap replicats of the gene tree (`qdist`, given for both 16S and 23S)
  
    ~~~bash
    ./TripVote/tq_distance_to_refs.py -i [bootstrap replicates] -r [main file] -m quartet -o [outputfile]
    ~~~~
  * the AAD range used to find this subset (`range`)
  Example lines look like this:

  |range|pdist|qdist16|qdist23|qdistAstral16|qdistAstral23|
  |-----|-----|-------|-------|-------------|-------------|
  |000_003|2.64952380952381|0.338385714285714|0.438738095238096|0.680952|0.514286|
  |000_003|2.53|0.409723809523809|0.266171428571429|0.52381|0.609524|

* [plot.R](plot.R): plots Figure 3c 
