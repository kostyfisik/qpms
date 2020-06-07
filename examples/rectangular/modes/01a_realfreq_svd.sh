#!/bin/bash
SCRIPTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
MISCDIR=../../../misc

source ${SCRIPTDIR}/00_params.sh

# try several lMaxes

${MISCDIR}/lat2d_realfreqsvd.py \
       	-B $BG_REFINDEX \
	-b s${A1X_nm}e-9 s${A1Y_nm}e-9 \
	-b s${A2X_nm}e-9 s${A2Y_nm}e-9 \
	-p s${P1X_nm}e-9 s${P1Y_nm}e-9 \
	-L 2 -m $METAL -r ${RADIUS_nm}e-9 -H ${HEIGHT_nm}e-9 \
	-k 0 0 \
	-F 2.001 0.001 2.250  \
	-P
