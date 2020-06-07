#!/bin/bash
SCRIPTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
MISCDIR=../../../misc

source ${SCRIPTDIR}/00_params.sh


for LMAX in 1 2 3 ; do  # try several cutoffs
  ${MISCDIR}/infiniterectlat-k0realfreqsvd.py \
       	-B $BG_REFINDEX \
	-p ${PX_nm}e-9 ${PY_nm}e-9 \
	-L $LMAX -m $METAL -r ${RADIUS_nm}e-9 -H ${HEIGHT_nm}e-9 \
	-F 2.001 0.001 2.250  \
	-P
done
