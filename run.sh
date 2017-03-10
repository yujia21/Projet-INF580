#!/bin/bash

for n in 10 20 30 40 50 60 70 80 90
do
    for p in 0.3 0.5 0.7 
    do
        DAT=data/${n}_${p}.edgelist
        python MaxCut.py $DAT erdos_renyi
    done
done

