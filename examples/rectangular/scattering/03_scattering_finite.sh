#!/bin/bash
SCRIPTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
MISCDIR=../../../misc

source ${SCRIPTDIR}/00_params.sh

for PSI in 1; do
  ${MISCDIR}/finiterectlat-scatter.py \
	--size 20 20 \
       	-B $BG_REFINDEX \
	-p ${PX_nm}e-9 ${PY_nm}e-9 \
	-L 3 -m $METAL -r ${RADIUS_nm}e-9 -H ${HEIGHT_nm}e-9 \
	--theta "s-0.005:0.005|101" \
	--phi 0 \
        --psi $PSI \
	--chi 0	\
	-P \
	-f "s2.150:2.180|100" \
	
done

