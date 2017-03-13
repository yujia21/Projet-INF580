#!/bin/bash

for k in 1 2 3 4
do
    for p in 0 1 2 3 4 
    do
        DAT=data_chimera.nonweighted.reduced2k2/${k}_${p}.edgelist
        python MaxCut.py $DAT chimera
    done
done

