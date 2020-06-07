#!/bin/bash
SCRIPTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
MISCDIR=../../../misc

source ${SCRIPTDIR}/00_params.sh

for bbb in 1 -2 -3 -4 ; do
  for coeff in $(seq 0.80 0.01 1.50 | sed -e s/,/./g) ; do
    ${MISCDIR}/lat2d_modes.py \
       	-n $BG_REFINDEX \
	-b s${A1X_nm}e-9 s${A1Y_nm}e-9 \
	-b s${A2X_nm}e-9 s${A2Y_nm}e-9 \
	-p s${P1X_nm}e-9 s${P1Y_nm}e-9 \
	-p s${P2X_nm}e-9 s${P2Y_nm}e-9 \
	-L 3 -m $METAL -r ${RADIUS_nm}e-9 -H ${HEIGHT_nm}e-9 \
	-k s$(bc <<< ${KPOINTX_nmi}*${coeff})e9 s$(bc <<< ${KPOINTY_nmi}*${coeff})e9 \
	-d $bbb \
	-t 1e13 \
	-T 0.2 \
	-c 250 
  done
done
