#!/bin/bash
#SBATCH --output=corr2.txt      # nom du fichier de sortie


for lh in 0.0000005 0.000001 0.000005 0.00001
do
    for lj in 0.00005 0.0001 0.0005 0.001
    do
        echo "lj $lj lh $lh"
        julia -t 16 test.jl --lambdaH $lh --lambdaJ $lj --mult 5 --rinit "../data/PF00014.fasta"
    done
done
