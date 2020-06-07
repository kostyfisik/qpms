#!/bin/bash

# Common parameters for a rectangular array
# N.B. We put those into bc, which does not understant exponential notation

export PX_nm=375
export PY_nm=375
export RADIUS_nm=30
export HEIGHT_nm=30
export METAL=Ag
export BG_REFINDEX=1.52


# Setup bc

echo 'scale=20;pi=3.14159265358979323846;' > bc_env
export BC_ENV_ARGS="bc_env"


# We have only one particle per unit cell here
export P1X_nm=0
export P1Y_nm=0


# Lattice vectors (for the general scripts)
export A1X_nm=${PX_nm}
export A1Y_nm=0
export A2X_nm=0
export A2Y_nm=${PY_nm}

# Reciprocal lattice vectors
export B1X_nmi=$(bc <<< '1/'${PX_nm})
export B1Y_nmi=0
export B2X_nmi=0
export B2Y_nmi=$(bc <<< '1/'${PY_nm})


