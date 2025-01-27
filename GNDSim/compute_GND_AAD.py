from sequence_lib import read_fasta, p_distance
from sys import argv

d=argv[1]
base_genome=argv[2]

# nt
names_mutated,seqs_mutated = read_fasta(d+"/mutated.fasta")
names_org,seqs_org = read_fasta(base_genome+"_contigs.fasta")

genome_org = {}
for name,seq in zip(names_org,seqs_org):
    genome_org[name] = seq

genome_mu = {}
for name,seq in zip(names_mutated,seqs_mutated):
    genome_mu[name] = seq

mutated = 0
for x in genome_org:
    m = p_distance(genome_org[x],genome_mu[x])
    # print(m)
    mutated += m*len(genome_org[x])

total = sum(len(s) for s in seqs_org)
gnd = mutated/total

#aa
names_mutated,seqs_mutated = read_fasta(d+"/mutated_genes.faa")
names_org,seqs_org = read_fasta(base_genome+"_contigs_genes.faa")

genome_org = {}
for name,seq in zip(names_org,seqs_org):
    genome_org[name] = seq

genome_mu = {}
for name,seq in zip(names_mutated,seqs_mutated):
    genome_mu[name] = seq

mutated = 0
for x in genome_org:
    m = p_distance(genome_org[x],genome_mu[x])
    # print(m)
    mutated += m*len(genome_org[x])

total = sum(len(s) for s in seqs_org)
aad = mutated/total

print(gnd,aad)
with open(d+"/gnd_aad.txt",'w') as fout:
    fout.write(str(gnd) + " " + str(aad))

