#!/bin/bash

DIR=data

n_list = [10,20,30,40,50,60,70,80,90]
p_list = [0.3,0.5,07]
for n in n_list : 
    for p in p_list : 
        DAT=$DIR/${n}_${p}.dat
        python MaxCut.py $DAT

