#!/bin/bash

## Assumes the availability of the newick utilities toolkit

for x in UPP.*_reduced.fasta.treefile; do 
 echo $x
 cat $x |sed -e "s/_[0-9]*//g"|nw_prune -v - `cat AAD-genomes.txt| tr '\n' ' '`|nw_distance -n -mm - > $x.dist
 for i in `seq 2 $(cat $x.dist|wc -l)`; do cut -f1,$i $x.dist |sort -k2n|head -n5|tr '\n' '\t'|cut -f 2,5,6,7,8,9,10; done |tee $x.dist.closest.txt
done
