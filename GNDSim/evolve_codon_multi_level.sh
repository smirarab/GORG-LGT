#! /bin/bash

# alpha=5.28
alpha=22
set -x
set -e

base_dir="genome_pairs/pair_0"
for g in {"AG-359-G18","AG-894-K15"}; do
for i in $(seq 1 2 70); do 
    echo $i $alpha
    python evolve.py `echo $i/100 | bc -l` $alpha "codon" $g $base_dir
    python compute_GND_AAD.py $base_dir/mutated_BLOSUM62_${g}_a$alpha\_aad$(printf "%03d" $i) $base_dir/$g
done
done

base_dir="genome_pairs/pair_1"
for g in {"AG-891-G05","AG-893-E23"}; do
for i in $(seq 1 2 70); do 
    echo $i $alpha
    python evolve.py `echo $i/100 | bc -l` $alpha "codon" $g $base_dir
    python compute_GND_AAD.py $base_dir/mutated_BLOSUM62_${g}_a$alpha\_aad$(printf "%03d" $i) $base_dir/$g
done
done

base_dir="genome_pairs/pair_2"
for g in {"AG-894-P05","AG-899-G06"}; do
for i in $(seq 1 2 70); do 
    echo $i $alpha
    python evolve.py `echo $i/100 | bc -l` $alpha "codon" $g $base_dir
    python compute_GND_AAD.py $base_dir/mutated_BLOSUM62_${g}_a$alpha\_aad$(printf "%03d" $i) $base_dir/$g
done
done

base_dir="genome_pairs/pair_3"
for g in {"AG-893-F23","AG-908-A02"}; do
for i in $(seq 1 2 70); do 
    echo $i $alpha
    python evolve.py `echo $i/100 | bc -l` $alpha "codon" $g $base_dir
    python compute_GND_AAD.py $base_dir/mutated_BLOSUM62_${g}_a$alpha\_aad$(printf "%03d" $i) $base_dir/$g
done
done

base_dir="genome_pairs/pair_4"
for g in {"AG-390-N04","AG-435-F03"}; do
for i in $(seq 1 2 70); do 
    echo $i $alpha
    python evolve.py `echo $i/100 | bc -l` $alpha "codon" $g $base_dir
    python compute_GND_AAD.py $base_dir/mutated_BLOSUM62_${g}_a$alpha\_aad$(printf "%03d" $i) $base_dir/$g
done
done
