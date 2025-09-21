#! /bin/bash

alpha=5

for g in {"AG-359-G18","AG-894-K15","AG-891-G05","AG-893-E23","AG-894-P05","AG-899-G06","AG-893-F23","AG-908-A02","AG-435-F03","AG-390-N04","AG-893-F11","AG-426-E17","AG-894-C14","AG-892-F15","AG-917-K06","AG-900-I13","AG-912-G22","AG-414-L04","AG-891-K05","AG-891-I18","AG-891-J07","AG-904-K03","AG-390-D15","AG-909-F14","AG-893-K09","AG-920-L07","AG-917-F18","AG-891-J04","AH-287-F15"}; do
echo $g
for i in {"0","5","10"} $(seq 20 10 90); do
    echo $i
    for s in {1..5}; do
    p=$(shuf -n 1 incompleteness.csv)
    ip=$(echo "$p * 100 / 1" | bc)
    python simulate_incompleteness.py -i simulated_genomes/$g/mutated_BLOSUM62_${g}_a${alpha}\_gnd$(printf "%05d" ${i})/mutated.fasta -p $p -s $s -o simulated_genomes/$g/mutated_BLOSUM62_${g}_a${alpha}\_gnd$(printf "%05d" ${i})/mutated_p${ip}_s${s}.fasta &
    done
done
for i in $(seq 1 1 20); do 
    echo $i
    for s in {1..5}; do
    p=$(shuf -n 1 incompleteness.csv)
    ip=$(echo "$p * 100 / 1" | bc)
    python simulate_incompleteness.py -i simulated_genomes/$g/mutated_BLOSUM62_${g}_a${alpha}\_gnd$(printf "%03d" ${i})/mutated.fasta -p $p -s $s -o simulated_genomes/$g/mutated_BLOSUM62_${g}_a${alpha}\_gnd$(printf "%03d" ${i})/mutated_p${ip}_s${s}.fasta &
    done
done
wait
done
