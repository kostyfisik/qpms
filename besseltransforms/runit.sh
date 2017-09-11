#!/bin/bash
K=$1
Q=$2
N=$3
module load mathematica
cat - vzor.m <<<"
kk=$K;
qq=$Q;
nn=$N;
" | math -noprompt > "${K}-${Q}-${N}"

