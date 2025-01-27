#! /usr/bin/env python
import statistics
from statistics import mean

from sequence_lib import read_fasta, write_fasta, p_distance
from random import *
from numpy.random import gamma, binomial, choice


def read_codon_table():
    codon = {}
    with open("codon_table.txt", 'r') as fin:
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
            for a2, p in zip(aa[1:], row[1:]):
                M[a1][a2] = float(p) if a2 != a1 else 0
    return M


def normalize(v):
    s = sum(v)
    return [x / s for x in v]


def randomize_rates(n_gene, alpha):
    rates = gamma(alpha, 1 / alpha, n_gene)
    # print("\nrates: ".join((str(r) for r in rates)))
    return normalize(rates)


def evolve(gseqs, g_weights, n_mus):
    # n_mus: number of mutations
    dist = {}
    print("Computing population ...")
    population = [(i, j) for i in range(len(gseqs)) for j in range(len(gseqs[i]))]
    print(len(population))

    print("Computing weights ...")
    weights = normalize([g_weights[i][j] for i in range(len(g_weights)) for j in range(len(g_weights[i]))])
    print("w:", sum(weights), n_mus, len(weights), len(set(weights)))

    print("Choosing mutation sites ...")
    mutated = choice(range(len(weights)), size=n_mus, replace=False, p=weights)
    geneDs = [[] for i in range(len(gseqs))]
    for k in mutated:
        i, j = population[k]
        geneDs[i].append(j)
        seq = (gseqs[i] + '.')[:-1]
        codon = (seq[(j // 3) * 3:(j // 3) * 3 + 3] + '.')[:-1]
        new_codon = "TAA"
        while new_codon in set(['TAA', 'TAG', 'TGA', 'ATG']):
            c = choice(list(set(['A', 'C', 'G', 'T']) - set([seq[j]])))
            new_codon = codon[0:j % 3] + c + codon[j % 3 + 1:3]
        gseqs[i] = seq[:j] + c + seq[j + 1:]
    # print("\nM: ".join(str(len(i)) for i in [[]]+geneDs))
    return geneDs


def evolve2(seqs, g_weights, n_mus):
    # n_mus: number of mutations
    dist = {}

    print("Computing weights ...")
    # weights = normalize([g_weights[i][j] for i in range(len(g_weights)) for j in range(len(g_weights[i]))])
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
        seqs[i] = seq[:j] + \
                  choice(list(set(['A', 'C', 'G', 'T']) - set([seq[j]]))) + \
                  seq[j + 1:]


def extract_locations(seqs, names, gnames):
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


def assign_weights(gseqs, g_locations, alpha):
    # randomly assign a (relative) mutation rate for each gene
    # the rate multipliers are drawn from a gamma distribution
    # all the loci in each gene are assigned the rate multipliers
    # of that gene as their weights; instead for the start and
    # end codons and the overlapping regions between the genes
    # are assigned weights 0 (i.e. not allowed to mutate)
    # alpha control the variance of the rate (i.e. shape of the
    # gamma distribution)
    n_genes = len(gseqs)
    g_rates = randomize_rates(n_genes, alpha)
    g_weights = []
    prv_node, prv_start, prv_end = (None, None, None)

    for rate, (node, start, end, rev) in zip(g_rates, g_locations):
        w = [0, 0, 0] + [rate] * (end - start + 1 - 6) + [0, 0, 0]
        # weight 0 for start and end codons
        g_weights.append(w)
        # assert len(w) == len(seq)
        # print(sum(g==0 for g in w))
    # print(sum(sum(pw for pw in w) for w in g_weights) / sum(len(w) for w in g_weights))
    return g_weights, [g / mean(g_rates) for g in g_rates]


def assign_weights2(gseqs, g_locations, alpha, seqs, names):
    # randomly assign a (relative) mutation rate for each gene
    # the rate multipliers are drawn from a gamma distribution
    # all the loci in each gene are assigned the rate multipliers
    # of that gene as their weights; instead for the start and
    # end codons and the overlapping regions between the genes
    # are assigned weights 0 (i.e. not allowed to mutate)
    # alpha control the variance of the rate (i.e. shape of the
    # gamma distribution)
    n_genes = len(gseqs)
    g_rates = randomize_rates(n_genes, alpha)
    g_weights = []
    prv_node, prv_start, prv_end = (None, None, None)

    for s in seqs:
        g_weights.append([None] * len(s))

    c = 0
    for rate, (node, start, end, rev) in zip(g_rates, g_locations):
        for j in [i for i in range(start, start + 3)] + [i for i in range(end - 2, end + 1)]:
            g_weights[names.index(node)][j] = 0
        for j in range(start + 3, end - 2):
            if g_weights[names.index(node)][j] == 0:
                continue
            if g_weights[names.index(node)][j] is None:
                g_weights[names.index(node)][j] = rate
            else:
                g_weights[names.index(node)][j] = choice([rate, g_weights[names.index(node)][j]], size=1)[0]
        # assert len(w) == len(seq)
        # print(sum(g==0 for g in w))
    # print(sum(sum(pw for pw in w) for w in g_weights) / sum(len(w) for w in g_weights))
    # print(g_weights)
    return g_weights, [g / mean(g_rates) for g in g_rates]


def reverse_complement(seq):
    # compute the DNA reverse complement
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    rseq = ''
    for x in seq[::-1]:
        rseq += complement[x]
    return rseq


codons = read_codon_table()
backcodon = dict()
for p, ns in codons.items():
    for n in ns:
        backcodon[n] = p


def translate(new_nt, rev):
    if rev:
        new_nt = reverse_complement(new_nt)
    prot = []
    for i in range(len(new_nt) // 3):
        new_codon = new_nt[3 * i:3 * i + 3]
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
        newAA = translate(genome[node][start:end + 1], rev)
        for i, aa in enumerate(newAA):
            if aa == "*" or i == 0 or i == len(newAA) - 1:
                if rev:
                    # print(start,end,end - i*3 - 3 + 1 ,end - i*3 + 1)
                    genome[node] = genome[node][0:end - i * 3 - 3 + 1] + \
                                   origs[names.index(node)][end - i * 3 - 3 + 1:end - i * 3 + 1] + \
                                   genome[node][end - i * 3 + 1:len(genome[node])]
                else:
                    genome[node] = genome[node][0:start + i * 3] + \
                                   origs[names.index(node)][start + i * 3:start + i * 3 + 3] + \
                                   genome[node][start + i * 3 + 3:len(genome[node])]
        new_aas.append(translate(genome[node][start:end + 1], rev))
    return list(genome.keys()), list(genome.values()), new_aas


def stitch_back(names, seqs, g_locations, gseqs, g_weights, gnames, Ds):
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
    for gseq, (node, start, end, rev), g_weight, gname, eD in zip(gseqs, g_locations, g_weights, gnames, Ds):
        ref = genome[node][start:end + 1]
        dna = reverse_complement(gseq) if rev else gseq
        new_aas.append(translate(dna, rev))
        D = sum(x != y for x, y in zip(dna, ref))
        if D != len(eD):
            print("*******", D, len(eD), node)
            print(" -- ".join("%d %s,%s" % (i, x, y) for i, (x, y) in enumerate(zip(dna, ref)) if x != y))
            print("     ", ", ".join(str(e) for e in sorted(eD)))
            # if D<len(eD):
            #    raise Exception("%d %d %s" %(D, eD, gname))
        # print("D:", D)
        if D / len(dna) <= 0.01:
            se += 1
        if 500 < end - start < 1500:
            # print("gene:r", D, g_weight)
            if D / len(dna) <= 0.01:
                ser += 1
            serc += 1
        ts += D
        gl += len(dna)
        new_seq = genome[node][:start] + dna + genome[node][end + 1:]
        genome[node] = new_seq
    print("sum:", se / len(gseqs), ser / serc, se, len(gseqs), ts, gl, ts / gl)
    return list(genome.keys()), list(genome.values()), new_aas


from sys import argv
from os import mkdir, getcwd, rmdir, listdir

p = float(argv[1])  # p is the proportion of aa mutations (i.e. 1-AAI)
alpha = float(argv[2])

names, seqs = read_fasta("AG-892-F15_contigs.fasta")  # full genome
gnames, _ = read_fasta("AG-892-F15_contigs_genes.faa")  # genes

g_locations, gseqs = extract_locations(seqs, names, gnames)
g_weights, g_rates = assign_weights2(gseqs, g_locations, alpha, seqs, names)
# g_weights, g_rates = assign_weights(gseqs, g_locations, alpha)

ls = sorted([end-start+1 for (node, start, end, rev) in g_locations])
print("Lengths", mean(ls), ls[len(g_locations)//2], statistics.variance(ls))

n_sites = sum(len(s) for s in gseqs)
n_mus = round(n_sites * p)
print(n_sites, n_mus)

# geneDs = evolve2(gseqs, g_weights, n_mus)
orig = [str(s) for s in seqs]
evolve2(seqs, g_weights, n_mus)
mutated_names, mutated_seqs, new_aaseqs = translate_all(names, seqs, g_locations, gnames, orig)


outdir = "mutated_a" + str(int(alpha)) + "_aad" + str(round(p * 100)).rjust(5, '0')
mkdir(outdir)
write_fasta(outdir + "/mutated_genes.faa", gnames, new_aaseqs)
write_fasta(outdir + "/mutated.fasta", mutated_names, mutated_seqs)
