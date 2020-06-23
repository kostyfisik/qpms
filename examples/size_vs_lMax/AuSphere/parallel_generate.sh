#!/bin/bash
#SBATCH --mem=200
#SBATCH -t 30:00
#SBATCH -c 4
#SBATCH -p short-ivb
#SBATCH --array=0-250

cat $0

contour_points=410

#radii_nm=(`seq 80 1 150`)
radii_nm=(`seq 50 1 300`)
radius_nm=${radii_nm[$SLURM_ARRAY_TASK_ID]}

for lMax in $(seq 1 5); do
  srun rectlat_simple_modes.py -p 580e-9 -m '4+0.7j' -r ${radius_nm}e-9 -k 0 0 --kpi -n 1.52 -L lMax -t 1e11 -b -2 -f 0.1 -i 1. -T .3 -N ${contour_points} 
done

