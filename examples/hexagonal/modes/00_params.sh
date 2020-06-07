#!/bin/bash
echo 'scale=20;pi=3.14159265358979323846;' > bc_env
export BC_ENV_ARGS="bc_env"

# We put those into bc, which does not understant exponential notation
SEPARATION_nm=576

# Particle positions within unit cell
export P1X_nm=0
export P1Y_nm=$(bc <<< ${SEPARATION_nm}/2)
export P2X_nm=0
export P2Y_nm=-$P1Y_nm

# Lattice vectors
export A1X_nm=$(bc <<< ${SEPARATION_nm}'*sqrt(3)')
export A1Y_nm=0
export A2X_nm=$(bc <<< ${SEPARATION_nm}'*sqrt(3)/2')
export A2Y_nm=$(bc <<< ${SEPARATION_nm}'*3/2')

# Reciprocal lattice vectors
export B1X_nmi=$(bc <<< '2*pi/sqrt(3)/'${SEPARATION_nm})
export B1Y_nmi=$(bc <<< '-2*pi/3/'${SEPARATION_nm})
export B2X_nmi=0
export B2Y_nmi=$(bc <<< '4*pi/3/'${SEPARATION_nm})

# a K-point coordinates
export KPOINTX_nmi=$(bc <<< '4*pi/3/sqrt(3)'/${SEPARATION_nm})
export KPOINTY_nmi=0.0 #$(bc <<< '4*pi/3/sqrt(3)'/${SEPARATION_nm})

# a M-point coordinates
export MPOINTX_nmi=0.0 
export MPOINTY_nmi=$(bc <<< '2*pi/3'/${SEPARATION_nm})


export RADIUS_nm=50
export HEIGHT_nm=50
export METAL=Au
export BG_REFINDEX=1.52
