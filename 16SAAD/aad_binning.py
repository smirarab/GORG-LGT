#! /usr/bin/env python
# Do the binning of the ANI distances based on the aad_tree constructed by TreeN93

from treeswift import *
from random import sample,choice,shuffle
import numpy as np

def compute_distance_from_leaves(tree):
    for node in tree.traverse_postorder():
        if node.is_leaf():
            node.t = 0
        else:
            c = node.child_nodes()[0]
            node.t = c.t + c.get_edge_length()

def two_way_cluster(tree,t1,t2):
# assumming tree is ultrametric
    C = []
    for node in tree.traverse_postorder():
        node.is_sub_cluster = False   
        for c in node.child_nodes():
            if node.t >= t2 and c.t < t2:
                C.append(c)
            if node.t > t1 and c.t <= t1:
                c.is_sub_cluster = True
    two_way_C = []
    for croot in C:
        C1 = []
        for node in croot.traverse_preorder():
            if node.is_sub_cluster:
                C1.append([ x.label for x in  node.traverse_leaves() ])         
        two_way_C.append(C1)                       
    return two_way_C

def list_pdist(D,taxalist,missing=60):
    pdist = []
    for i,x in enumerate(taxalist):
        for y in taxalist[i+1:]:
            key = (min(x,y),max(x,y))
            pdist.append(D[key] if key in D else missing)
    return pdist 


aad_tree = read_tree_newick("aad_tree_pruned.nwk")
compute_distance_from_leaves(aad_tree)

t2 = 1
t1 = 0
prev_t1 = -1
prev_t2 = 0
n = 10
k = 50 
win_step = 1
search_step = 1

binning = {}

while t2 <= 60:
    two_way_C = two_way_cluster(aad_tree,t1,t2)
    valid_C = []
    taxaList = []
    count = 0
    for C in two_way_C:
        if len(C) >= n:
            valid_C.append(C)
    if valid_C:
        while count < k:
            C = choice(valid_C)         
            c_indices = sample(list(range(len(C))),n)   
            L = []
            for idx in c_indices:
                L.append(choice(C[idx]))
            taxaList.append(L)
            count += 1
        binning[(t1,t2)] = taxaList
        prev_t1 = t1
        prev_t2 = t2
        t1 = t1 + win_step
        t2 = t2 + search_step
    elif t2 >= 60:
        binning.pop((prev_t1,prev_t2))    
        t1 = prev_t1
    else:
        t2 += search_step

tree_16S = read_tree_newick("RAxML.T16S-SINA-1000-bootcutoff_MVrooted_pruned")
tree_23S = read_tree_newick("RAxML.T23S-SINA-1000-bootcutoff_MVrooted_pruned")
tree_astral = read_tree_newick("astral_pruned.tre")

D = {}
with open("aad_distances.txt",'r') as fin:
    for line in fin:
        s1,s2,d = line.strip().split(',')
        D[(min(s1,s2),max(s1,s2))] = float(d)

for (t1,t2) in binning:
    file_pdist = "pdist_" + str(t1).rjust(3,'0') + "_" + str(t2).rjust(3,'0') + ".txt"
    file_16S = "16S_" + str(t1).rjust(3,'0') + "_" + str(t2).rjust(3,'0') + ".nwk"
    file_23S = "23S_" + str(t1).rjust(3,'0') + "_" + str(t2).rjust(3,'0') + ".nwk"
    file_astral = "astral_" + str(t1).rjust(3,'0') + "_" + str(t2).rjust(3,'0') + ".nwk"
    fout_pdist = open(file_pdist,'w')
    fout_16S = open(file_16S,'w')
    fout_23S = open(file_23S,'w')
    fout_astral = open(file_astral,'w')

    pdist = []
    for taxaList in binning[(t1,t2)]:
        pdist += list_pdist(D,taxaList,missing=60)
        fout_16S.write(tree_16S.extract_tree_with(taxaList).newick()+"\n")
        fout_23S.write(tree_23S.extract_tree_with(taxaList).newick()+"\n")
        fout_astral.write(tree_astral.extract_tree_with(taxaList).newick()+"\n")
   
    for d in pdist:
        fout_pdist.write(str(d)+"\n")
    
    fout_16S.close()    
    fout_23S.close()    
    fout_astral.close()    
    fout_pdist.close()
