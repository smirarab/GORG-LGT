import sys
import os
import argparse
import time
from sequence_lib import read_fasta, write_fasta, p_distance
import numpy as np

def add_missing(names, seqs, p, C):
	L = sum([len(s) for s in seqs])
	num_regions = np.random.poisson(lam = C)
	for _ in range(num_regions):
		length = np.random.poisson(lam = p*L/num_regions)
		if  length == 0:\
			continue
		seq_pos = [(i, j) for i in range(len(seqs)) for j in range(len(seqs[i]) - length + 1)]

		pos = seq_pos[np.random.choice(len(seq_pos))]

		s = seqs[pos[0]]
		s1 = s[:pos[1]]
		s2 = s[pos[1] + length:]

		# print(len(s), length, len(s1), len(s2))

		seqs[pos[0]] = s1

		if len(s2) == 0:
			continue
			
		seqs.insert(pos[0] + 1, s2)

		n = names[pos[0]]

		new_id = int(n.split("_")[-1]) + 1
		new_n = "_".join(n.split("_")[:-1]) + "_" + str(new_id)

		while new_n in names:
			new_id += 1
			new_n = "_".join(n.split("_")[:-1]) + "_" + str(new_id)

		names.insert(pos[0] + 1, new_n)

		# print(names[pos[0]], len(seqs[pos[0]]))
		# print(names[pos[0]+1], len(seqs[pos[0]+1]))
	# print(*names, sep = "\n")
	return names, seqs



def main():
	parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-i', '--input', required=True, help="Input fasta file")
	parser.add_argument('-p', '--incompleteness', required=True, help="Incompleteness")
	parser.add_argument('-s', '--seed', required=True, help="Random Seed")
	parser.add_argument('-o', '--output', required=True, help="Output file")

	args = parser.parse_args()

	np.random.seed(int(args.seed))


	p = float(args.incompleteness)
	C = 17

	names,seqs = read_fasta(args.input)
	# print(len(names))

	new_names, new_seqs = add_missing(names, seqs, p, C)
	# print(len(new_names))
	write_fasta(args.output, new_names, new_seqs)
	# print(seqs)

if __name__ == "__main__":
	main()