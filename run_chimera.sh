#!/bin/bash

for k in 5 6 7 8
do
    for p in 0
    do
        DAT=data_chimera.nonweighted/${k}_${p}.edgelist
        python MaxCutVis.py $DAT chimera
    done
done

