Code and data to reproduce Several figures in the paper related to gene trees. 

## Requirement:

* The bash script `comptue-closest.sh` relies on the fantastic [newick utilities](https://github.com/tjunier/newick_utils) toolkit. But results of it are provided.
* R Scripts require R (works with 4.2.2)

## Files

* [GORGv1_16SSAGs_aai_summary.csv.xz](GORGv1_16SSAGs_aai_summary.csv.xz): Pairwise AAD results

* [all-gene-trees.zip](all-gene-trees.zip): all gene family trees in newick format

* [comptue-closest.sh](comptue-closest.sh): a script to compute the closest distance from the gene trees. It relies on the fantastic [newick utilities](https://github.com/tjunier/newick_utils) toolkit.

* `all-closest-dist.txt.xz`: top three closest leaves to each leaf in the gene tree (query) and their path length on the gene tree to the query

  |gene|queryGenome|closest1|closest1BL|closest2|closest2BL|closest3|closest3AAD|
  |-|-|-|-|-|-|-|-|
  |accA|AG-359-D09|AH-321-P04|0.179366|AH-321-K15|0.455919|AG-337-G21|0.484358|
  |accA|AH-321-P04|AG-359-D09|0.179366|AH-321-K15|0.472453|AG-337-G21|0.500892|

  This file is created using the following commands.
   
  ~~~bash
  unzip all-gene-trees.zip
  bash comptue-closest.sh
  grep "" UPP.*_reduced.fasta.treefile.dist.closest.txt|sed -e "s/.*UPP.//" -e "s/_reduced.fasta.treef.*:/\t/g" > all-closest-dist.txt
  ~~~
* `minimum-AAD-per-gene.csv.xz`: For each gene, it includes the minimum possible AAD to each leaf present in that gene tree *among all the genomes present in that gene tree*

* [plot.R](plot.R): code to plot the data above (Figures 2B, S1, S4, and S7).

* [outlier.R](plot.R): code to reproduce [hellinger_results.csv](hellinger_results.csv).

* [multiplicities.txt](multiplicities.txt): The multiplicity of species per gene family tree 

* `AAD-genomes.txt`: List of genomes with AAD values computed and used
