#!/bin/bash
SCRIPTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
MISCDIR=../../../misc

source ${SCRIPTDIR}/00_params.sh

for PSI in 0 1; do
  ${MISCDIR}/infiniterectlat-scatter.py \
       	-B $BG_REFINDEX \
	-p ${PX_nm}e-9 ${PY_nm}e-9 \
	-L 3 -m $METAL -r ${RADIUS_nm}e-9 -H ${HEIGHT_nm}e-9 \
	--theta "s-0.015:0.015|201" \
	--phi 0 \
        --psi $PSI \
	--chi 0	\
	-f "s2.110:2.230|100" \
	-P
done

