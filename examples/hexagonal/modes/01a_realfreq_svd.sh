#!/bin/bash
SCRIPTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
MISCDIR=../../../misc

source ${SCRIPTDIR}/00_params.sh

${MISCDIR}/lat2d_realfreqsvd.py \
       	-B $BG_REFINDEX \
	-b s${A1X_nm}e-9 s${A1Y_nm}e-9 \
	-b s${A2X_nm}e-9 s${A2Y_nm}e-9 \
	-p s${P1X_nm}e-9 s${P1Y_nm}e-9 \
	-p s${P2X_nm}e-9 s${P2Y_nm}e-9 \
	-L 3 -m $METAL -r ${RADIUS_nm}e-9 -H ${HEIGHT_nm}e-9 \
	-k s${KPOINTX_nmi}e9 s${KPOINTY_nmi}e9 \
	-F 1.3 0.001 1.5 \
	-P
