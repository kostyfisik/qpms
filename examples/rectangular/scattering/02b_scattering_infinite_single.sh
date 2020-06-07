#!/bin/bash
SCRIPTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
MISCDIR=../../../misc

source ${SCRIPTDIR}/00_params.sh

for PSI in 0 1; do
  ${MISCDIR}/infiniterectlat-scatter.py \
       	-B $BG_REFINDEX \
	-p ${PX_nm}e-9 ${PY_nm}e-9 \
	-L 2 -m $METAL -r ${RADIUS_nm}e-9 -H ${HEIGHT_nm}e-9 \
	--theta "s-0.05:0.05|101" \
	--phi 0 \
        --psi $PSI \
	--chi 0	\
	-P \
	-f 2.15
done

