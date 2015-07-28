#!/usr/bin/env bash

./yeast_unequal_ascus 0 1 40 5 0.1 0.8
#./yeast 0.05 0

flist=(`ls sim*`)
length=${#flist[@]}
lastposition=$((length-1))
./plot_sims.r ${flist[${lastposition}]}

