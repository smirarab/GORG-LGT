#! /usr/bin/env python

from sequence_lib import read_fasta, write_fasta
from random import *
from numpy.random import gamma, binomial, choice
import numpy as np
import json

def read_codon_table():
    codon = {}
    with open("codon_table.txt",'r') as fin:
        for line in fin:
            code = line.strip().split()
            codon[code[0]] = code[1:]
    return codon

def read_blosum62():
    M = {}
    with open("blosum62sym.csv") as fin:
        aa = [x[1] for x in fin.readline().strip().split(",")]
        for line in fin:
            row = line.strip().split(",")
            a1 = row[0][1]
            M[a1] = {}
            for a2,p in zip(aa[1:],row[1:]):
                M[a1][a2] = float(p) if a2 != a1 else 0
    return M

def normalize(v):
    s = sum(v)
    return [x/s for x in v]

def randomize_rates(n_gene,alpha):    
    rates = gamma(alpha,1/alpha,n_gene)
    return normalize(rates)

def evolve(gseqs,g_weights,n_mus,M):
    # n_mus: number of mutations

    print("Computing population ...")
    population = [(i,j) for i in range(len(gseqs)) for j in range(len(gseqs[i]))]

    print("Computing weights ...")
    weights = normalize([g_weights[i][j] for i in range(len(g_weights)) for j in range(len(g_weights[i]))])

    print("Choosing mutation sites ...")
    mutated = choice(range(len(weights)),size=n_mus,replace=False,p=weights)

    for k in mutated:
        i,j = population[k]
        seq = gseqs[i]
        c = choice(list(M[seq[j]].keys()),p=normalize(M[seq[j]].values()))
        gseqs[i] = seq[:j] + c + seq[j+1:]

def evolve2(seqs, g_weights, n_mus):
    # n_mus: number of mutations
    dist = {}

    print("Computing weights ...")
    #weights = normalize([g_weights[i][j] for i in range(len(g_weights)) for j in range(len(g_weights[i]))])
    weights = normalize([g if g is not None else 0 for gs in g_weights for g in gs])
    print("w:", sum(weights), n_mus, len(weights), len(set(weights)))

    print("Choosing mutation sites ...")
    mutated = choice(range(len(weights)), size=n_mus, replace=False, p=weights)
    i = 0
    rs = len(seqs[i])
    for k in sorted(mutated):
        while k > rs:
            i += 1
            rs += len(seqs[i])
        seq = seqs[i]
        j = k - rs + len(seqs[i])
        seqs[i] = seq[:j] +\
                  choice(list(set(['A', 'C', 'G', 'T']) - set([seq[j]]))) + \
                  seq[j + 1:]

def extract_locations(gnames):
# parse the gene names to extract the start and end points
    g_locations = []    
    #g_locations_adjusted = []
    #i_locations = []

    for i,g in enumerate(gnames):
        node,start,end,direction = g.strip().split("#")[:4]
        node = "_".join(node.strip().split("_")[:-1])
        start = int(start)-1
        end = int(end)-1
        rev = float(direction) < 0
        g_locations.append((node,start,end,rev))
    
    return g_locations

def extract_locations2(seqs, names, gnames):
    # parse the gene names to extract the start and end points
    g_locations = []
    # g_locations_adjusted = []
    # i_locations = []
    gseqs = []
    ends = set()
    for i, g in enumerate(gnames):
        node, start, end, direction = g.strip().split("#")[:4]
        node = "_".join(node.strip().split("_")[:-1])
        start = int(start) - 1
        end = int(end) - 1
        rev = float(direction) < 0
        g_locations.append((node, start, end, rev))
        if rev:
            gseqs.append(reverse_complement(seqs[names.index(node)][start:end + 1]))
        else:
            gseqs.append(seqs[names.index(node)][start:end + 1])

    return g_locations, gseqs

def assign_weights(gseqs,g_locations,alpha):
# randomly assign a (relative) mutation rate for each gene
# the rate multipliers are drawn from a gamma distribution
# all the loci in each gene are assigned the rate multipliers 
# of that gene as their weights; instead for the start and 
# end codons and the overlapping regions between the genes 
# are assigned weights 0 (i.e. not allowed to mutate)
# alpha control the variance of the rate (i.e. shape of the
# gamma distribution) 
    n_genes = len(gseqs)
    # print(n_genes)
    g_rates = randomize_rates(n_genes,alpha)
    g_weights = []
    prv_node, prv_start, prv_end = (None,None,None)

    for seq,rate,(node,start,end,_) in zip(gseqs,g_rates,g_locations):
        w = [rate]*len(seq)
        # weight 0 for start and end codons
        w[0] = 0 
        w[-1] = 0
        # check for overlapping with previous genes
        # here we assume that the genes are ordered
        # and a gene can only overlap with at most
        # one other gene
        # print(prv_node, prv_start, prv_end)
        if prv_node == node:
            if prv_end > start:
                prv_w = g_weights[-1]
                # set zero-weights for prv_w
                e = prv_end
                i = -1
                while e > start:
                    prv_w[i] = 0
                    e -= 3
                    i -= 1
                # set zero-weights for w    
                s = start
                i = 1
                while s < prv_end:
                    w[i] = 0
                    s += 3
                    i += 1                
        # add w to g_weights
        g_weights.append(w) 

    # w = np.array(g_weights[0])
    # print(w[w == 0])

    all_pos = [len(g) for g in g_weights]
    zero_pos = np.array([len(np.array(g)[np.array(g) == 0]) for g in g_weights])
    # print(zero_pos[zero_pos != 2])

    # print(sum(zero_pos), sum(all_pos))
    # print(len(g_weights))  

    return g_weights

def assign_weights2(gseqs, g_locations, alpha, seqs, names, wfile):
    # randomly assign a (relative) mutation rate for each gene
    # the rate multipliers are drawn from a gamma distribution
    # all the loci in each gene are assigned the rate multipliers
    # of that gene as their weights; instead for the start and
    # end codons and the overlapping regions between the genes
    # are assigned weights 0 (i.e. not allowed to mutate)
    # alpha control the variance of the rate (i.e. shape of the
    # gamma distribution)
    if os.path.isfile(wfile):
        with open(wfile, 'r') as f:
            g_weights = json.load(f)
            print("loaded")
            return g_weights
    n_genes = len(gseqs)
    g_rates = randomize_rates(n_genes, alpha)
    g_weights = []
    prv_node, prv_start, prv_end = (None, None, None)

    for s in seqs:
        g_weights.append([None] * len(s))

    c = 0
    for rate, (node, start, end, rev) in zip(g_rates, g_locations):
        for j in [i for i in range(start,start+3)] + [i for i in range(end - 2, end +1)]:
            g_weights[names.index(node)][j] = 0
        for j in range(start + 3, end - 2):
            if g_weights[names.index(node)][j] == 0:
                continue
            if g_weights[names.index(node)][j] is None:
                g_weights[names.index(node)][j] = rate
            else:
                g_weights[names.index(node)][j] = choice([rate, g_weights[names.index(node)][j]], size=1)[0]
        #assert len(w) == len(seq)
        # print(sum(g==0 for g in w))
    #print(sum(sum(pw for pw in w) for w in g_weights) / sum(len(w) for w in g_weights))
    #print(g_weights)
    # print(g_weights)
    with open(wfile, 'w') as f:
        json.dump(g_weights, f)
    return g_weights
    # return g_weights, [g / mean(g_rates) for g in g_rates]

def reverse_complement(seq):
# compute the DNA reverse complement
    complement = {'A':'T','T':'A','G':'C','C':'G'}
    rseq = ''
    for x in seq[::-1]:
        rseq += complement[x]
    return rseq

def back_translate(peptide,ref,rev):
# back translation is not unique
# here we do the following:
# if a codon exactly matches the original, 
# then use it. Otherwise, randomly choose a codon 
# with the weight is the probability 
# that the reference dna mutated into that codon 
# eg: prob(ACT --> ACG) = 1/4*(3/4)**2; prob(ACT --> AAG) = (1/4)**2*(3/4)
    if rev:
        ref = reverse_complement(ref)
    codon = read_codon_table()
    dna = ''
    for i,aa in enumerate(peptide):
        ref_i = ref[3*i:3*i+3]
        w_total = 0
        translated = None
        # print("="*200)
        # print(ref_i)
        for tr in codon[aa]:
            d = sum(a != b for a,b in zip(tr,ref_i))
            if d == 0:
                translated = tr
                break
            #w = (0.25**d)*(0.75**(3-d))
        if not translated:
            for tr in codon[aa]:
                d = sum(a != b for a,b in zip(tr,ref_i))
                w = 4-d
                w_total += w
                if random() <= w/w_total: # choose this codon with probability w/w_total
                    translated = tr
                    break
        if translated != ref_i:
            pass
            # print(translated, ref_i)
            # print(sum(a != b for a,b in zip(translated,ref_i)))
        # print("=" * 200)
        dna += translated
    D = sum(x!=y for x,y in zip(dna,ref))/len(dna)
    if rev:
        dna = reverse_complement(dna)
    return dna,D 

def stitch_back(names,seqs,g_locations,gseqs):
    # hash genome names
    genome = {}
    for node,seq in zip(names,seqs):
        genome[node] = seq

    # stitch the genes to the genome
    for gseq,(node,start,end,rev) in zip(gseqs,g_locations):
        ref = genome[node][start:end+1]
        dna,D = back_translate(gseq,genome[node][start:end+1],rev)
        new_seq = genome[node][:start] + dna + genome[node][end+1:]
        genome[node] = new_seq

    return list(genome.keys()),list(genome.values())


def translate(new_nt, rev):
    if rev:
        new_nt = reverse_complement(new_nt)
    prot = []
    for i in range(len(new_nt)//3):
        new_codon = new_nt[3*i:3*i+3]
        prot.append(backcodon[new_codon])
    prot = "".join(prot)
    return prot

def translate_all(names, seqs, g_locations, gnames, origs):
    # hash genome names
    genome = {}
    new_aas = []
    for node, seq in zip(names, seqs):
        genome[node] = seq

    # stitch the genes to the genome
    se = 0
    ser = 0
    serc = 0
    ts = 0
    gl = 0
    for (node, start, end, rev), gname in zip(g_locations, gnames):
        newAA = translate(genome[node][start:end + 1],rev)
        for i, aa in enumerate(newAA):
            if aa == "*" or i == 0 or i == len(newAA)-1:
                if rev:
                    #print(start,end,end - i*3 - 3 + 1 ,end - i*3 + 1)
                    genome[node] = genome[node][0:end - i*3 - 3 + 1] + \
                                   origs[names.index(node)][end - i*3 - 3 + 1:end - i*3 +1] + \
                                   genome[node][end - i*3 + 1:len(genome[node])]
                else:
                    genome[node] =             genome[node][0:start+i*3] + \
                                   origs[names.index(node)][start+i*3:start+i*3+3] + \
                                               genome[node][start+i*3+3:len(genome[node])]
        new_aas.append(translate(genome[node][start:end + 1],rev))
    return list(genome.keys()), list(genome.values()), new_aas

def compute_n_mus(n_sites,p,p_g=0.95):
# n_sites is the number of DNA sites, while
# n_mus is the number of aa mutations
# p is the proportion of DNA mutations in the genome
# p_g is the proportion of genes to the full genome
    nt2aa = 1/2
    return round(n_sites*p*nt2aa)

from sys import argv
import os
from os import mkdir,getcwd,rmdir,listdir

p = float(argv[1]) # p is the proportion of aa mutations (i.e. 1-AAI)
alpha = float(argv[2])
# pair = argv[3]
model = argv[3]
base_genome = argv[4]
if model == "codon":
    base_dir = argv[5]
elif model == "nt":
    base_dir = "simulated_genomes/" + base_genome

M = read_blosum62()
names,seqs = read_fasta(base_dir + "/" + base_genome +"_contigs.fasta") # full genome
gnames,gseqs = read_fasta(base_dir + "/" + base_genome +"_contigs_genes.faa") # genes

#n_sites = sum(len(s) for s in seqs) 
#n_mus = compute_n_mus(n_sites,p)
n_sites = sum(len(s) for s in gseqs) 
n_mus = round(n_sites*p)
print(p, n_mus)

if model == 'codon':
    g_locations = extract_locations(gnames)
    g_weights = assign_weights(gseqs,g_locations,alpha)
    evolve(gseqs,g_weights,n_mus,M)
    mutated_names, mutated_seqs = stitch_back(names,seqs,g_locations,gseqs)

    if p < 1e-2:
        outdir = base_dir + "/mutated_BLOSUM62_" + base_genome + "_a" + str(int(alpha)) + "_aad" + str(round(p*10000)).rjust(5,'0')
    else:
        outdir = base_dir + "/mutated_BLOSUM62_" + base_genome + "_a" + str(int(alpha)) + "_aad" + str(round(p*100)).rjust(3,'0')
    if not os.path.isdir(outdir):
        mkdir(outdir)
    write_fasta(outdir+"/mutated_genes.faa",gnames,gseqs)
    write_fasta(outdir+"/mutated.fasta",mutated_names,mutated_seqs)

elif model == 'nt':
    weight_file = base_dir + "/" + "g_weights_a" + str(int(alpha)) + ".npy"

    g_locations, gseqs = extract_locations2(seqs, names, gnames)
    g_weights = assign_weights2(gseqs, g_locations, alpha, seqs, names, weight_file)
    #g_weights, g_rates = assign_weights(gseqs, g_locations, alpha)


    #geneDs = evolve2(gseqs, g_weights, n_mus)
    codons = read_codon_table()
    backcodon = dict()
    for ps,ns in codons.items():
        for n in ns:
            backcodon[n] = ps
    orig = [str(s) for s in seqs]
    evolve2(seqs, g_weights, n_mus)
    mutated_names, mutated_seqs, new_aaseqs = translate_all(names, seqs, g_locations, gnames, orig)

    if p < 1e-2:
        outdir = base_dir + "/mutated_BLOSUM62_" + base_genome + "_a" + str(int(alpha)) + "_gnd" + str(round(p*10000)).rjust(5,'0')
    else:
        outdir = base_dir + "/mutated_BLOSUM62_" + base_genome + "_a" + str(int(alpha)) + "_gnd" + str(round(p*100)).rjust(3,'0')
    if not os.path.isdir(outdir):
        mkdir(outdir)
    write_fasta(outdir + "/mutated_genes.faa", gnames, new_aaseqs)
    write_fasta(outdir + "/mutated.fasta", mutated_names, mutated_seqs)

