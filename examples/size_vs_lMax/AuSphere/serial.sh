#!/bin/bash

kx=0.0
contour_points=410

radii_nm=(`seq 50 1 300`)
radius_nm=${radii_nm[$SLURM_ARRAY_TASK_ID]}

for lMax in $(seq 1 5) ; do
  for radius_nm in $(seq 50 1 300) ; do
    rectlat_simple_modes.py -p 580e-9 -m 'Au' -r ${radius_nm}e-9 -k $kx 0 --kpi -n 1.52 -L $lMax -t 1e11 -b -2 -f 0.1 -i 1. -T .3 -N ${contour_points} --lMax-extend 10
  done 
done
