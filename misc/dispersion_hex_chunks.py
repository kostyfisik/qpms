#!/usr/bin/env python3

import argparse, re, random, string
import subprocess
from scipy.constants import hbar, e as eV, pi, c

def make_action_sharedlist(opname, listname):
    class opAction(argparse.Action):
        def __call__(self, parser, args, values, option_string=None):
            if (not hasattr(args, listname)) or getattr(args, listname) is None:
                setattr(args, listname, list())
            getattr(args,listname).append((opname, values))
    return opAction

parser = argparse.ArgumentParser()
#TODO? použít type=argparse.FileType('r') ?
parser.add_argument('--TMatrix', action='store', required=True, help='Path to TMatrix file')
#parser.add_argument('--griddir', action='store', required=True, help='Path to the directory with precalculated translation operators')
parser.add_argument('--output_prefix', action='store', required=True, help='Prefix to the npz output (will be appended frequency, hexside and chunkno)')
#sizepar = parser.add_mutually_exclusive_group(required=True)
parser.add_argument('--hexside', action='store', type=float, required=True, help='Lattice hexagon size length')
parser.add_argument('--plot_TMatrix', action='store_true', help='Visualise TMatrix on the first page of the output')
#parser.add_argument('--SVD_output', action='store', help='Path to output singular value decomposition result')
parser.add_argument('--maxlayer', action='store', type=int, default=100, help='How far to sum the lattice points to obtain the dispersion')
parser.add_argument('--scp_to', action='store', metavar='N', type=str, help='SCP the output files to a given destination')
parser.add_argument('--background_permittivity', action='store', type=float, default=1., help='Background medium relative permittivity (default 1)')
parser.add_argument('--eVfreq', action='store', required=True, type=float, help='Frequency in eV')
parser.add_argument('--kdensity', action='store', type=int, default=33, help='Number of k-points per x-axis segment')
parser.add_argument('--chunklen', action='store', type=int, default=1000, help='Number of k-points per output file (default 1000)')
parser.add_argument('--lMax', action='store', type=int, help='Override lMax from the TMatrix file')
#TODO some more sophisticated x axis definitions
parser.add_argument('--gaussian', action='store', type=float, metavar='σ', help='Use a gaussian envelope for weighting the interaction matrix contributions (depending on the distance), measured in unit cell lengths (?) FIxME).')
parser.add_argument('--verbose', '-v', action='count', help='Be verbose (about computation times, mostly)')
popgrp=parser.add_argument_group(title='Operations')
popgrp.add_argument('--tr', dest='ops', action=make_action_sharedlist('tr', 'ops'), default=list()) # the default value for dest can be set once
popgrp.add_argument('--tr0', dest='ops', action=make_action_sharedlist('tr0', 'ops'))
popgrp.add_argument('--tr1', dest='ops', action=make_action_sharedlist('tr1', 'ops'))
popgrp.add_argument('--sym', dest='ops', action=make_action_sharedlist('sym', 'ops'))
popgrp.add_argument('--sym0', dest='ops', action=make_action_sharedlist('sym0', 'ops'))
popgrp.add_argument('--sym1', dest='ops', action=make_action_sharedlist('sym1', 'ops'))
#popgrp.add_argument('--mult', dest='ops', nargs=3, metavar=('INCSPEC', 'SCATSPEC', 'MULTIPLIER'), action=make_action_sharedlist('mult', 'ops'))
#popgrp.add_argument('--mult0', dest='ops', nargs=3, metavar=('INCSPEC', 'SCATSPEC', 'MULTIPLIER'), action=make_action_sharedlist('mult0', 'ops'))
#popgrp.add_argument('--mult1', dest='ops', nargs=3, metavar=('INCSPEC', 'SCATSPEC', 'MULTIPLIER'), action=make_action_sharedlist('mult1', 'ops'))
popgrp.add_argument('--multl', dest='ops', nargs=3, metavar=('INCL[,INCL,...]', 'SCATL[,SCATL,...]', 'MULTIPLIER'), action=make_action_sharedlist('multl', 'ops'))
popgrp.add_argument('--multl0', dest='ops', nargs=3, metavar=('INCL[,INCL,...]', 'SCATL[,SCATL,...]', 'MULTIPLIER'), action=make_action_sharedlist('multl0', 'ops'))
popgrp.add_argument('--multl1', dest='ops', nargs=3, metavar=('INCL[,INCL,...]', 'SCATL[,SCATL,...]', 'MULTIPLIER'), action=make_action_sharedlist('multl1', 'ops'))
parser.add_argument('--frequency_multiplier', action='store', type=float, default=1., help='Multiplies the frequencies in the TMatrix file by a given factor.')
# TODO enable more flexible per-sublattice specification
pargs=parser.parse_args()
print(pargs)

maxlayer=pargs.maxlayer
hexside=pargs.hexside
eVfreq = pargs.eVfreq
freq = eVfreq*eV/hbar
verbose=pargs.verbose

TMatrix_file = pargs.TMatrix

epsilon_b = pargs.background_permittivity #2.3104
gaussianSigma = pargs.gaussian if pargs.gaussian else None # hexside * 222 / 7
interpfreqfactor = pargs.frequency_multiplier
scp_dest = pargs.scp_to if pargs.scp_to else None
kdensity = pargs.kdensity
chunklen = pargs.chunklen

ops = list()
opre = re.compile('(tr|sym|copy|multl|mult)(\d*)')
for oparg in pargs.ops:
    opm = opre.match(oparg[0])
    if opm:
        ops.append(((opm.group(2),) if opm.group(2) else (0,1), opm.group(1), oparg[1]))
    else:
        raise # should not happen
print(ops)


# -----------------finished basic CLI parsing (except for op arguments) ------------------
from qpms.timetrack import _time_b, _time_e
btime=_time_b(verbose)

import qpms
import numpy as np
import os, sys, warnings, math
from scipy import interpolate
nx = None
s3 = math.sqrt(3)


# specifikace T-matice zde
cdn = c/ math.sqrt(epsilon_b)
TMatrices_orig, freqs_orig, freqs_weirdunits_orig, lMaxTM = qpms.loadScuffTMatrices(TMatrix_file)
lMax = lMaxTM
if pargs.lMax:
    lMax = pargs.lMax if pargs.lMax else lMaxTM    
my, ny = qpms.get_mn_y(lMax)
nelem = len(my)
if pargs.lMax: #force commandline specified lMax
    TMatrices_orig = TMatrices_orig[...,0:nelem,:,0:nelem]

TMatrices = np.array(np.broadcast_to(TMatrices_orig[:,nx,:,:,:,:],(len(freqs_orig),2,2,nelem,2,nelem)) )

#TMatrices[:,:,:,:,:,ny==3] *= factor13inc
#TMatrices[:,:,:,ny==3,:,:] *= factor13scat
xfl = qpms.xflip_tyty(lMax)
yfl = qpms.yflip_tyty(lMax)
zfl = qpms.zflip_tyty(lMax)
c2rot = qpms.apply_matrix_left(qpms.yflip_yy(3),qpms.xflip_yy(3),-1)

reCN = re.compile('(\d*)C(\d+)')
#TODO C nekonečno

for op in ops:
    if op[0] == 'all':
        targets = (0,1)
    elif isinstance(op[0],int):
        targets = (op[0],)
    else:
        targets = op[0]
        
    if op[1] == 'sym':
        mCN = reCN.match(op[2]) # Fuck van Rossum for not having assignments inside expressions
        if op[2] == 'σ_z':
            for t in targets:
                TMatrices[:,t] = (TMatrices[:,t] + qpms.apply_ndmatrix_left(zfl,qpms.apply_ndmatrix_left(zfl, TMatrices[:,t], (-4,-3)),(-2,-1)))/2
        elif op[2] == 'σ_y':
            for t in targets:
                TMatrices[:,t] = (TMatrices[:,t] + qpms.apply_ndmatrix_left(yfl,qpms.apply_ndmatrix_left(yfl, TMatrices[:,t], (-4,-3)),(-2,-1)))/2
        elif op[2] == 'σ_x':
            for t in targets:
                TMatrices[:,t] = (TMatrices[:,t] + qpms.apply_ndmatrix_left(xfl,qpms.apply_ndmatrix_left(xfl, TMatrices[:,t], (-4,-3)),(-2,-1)))/2
        elif op[2] == 'C2': # special case of the latter
            for t in targets:
                TMatrices[:,t] = (TMatrices[:,t] + qpms.apply_matrix_left(c2rot,qpms.apply_matrix_left(c2rot, TMatrices[:,t], -3),-1))/2                
        elif mCN:
            rotN = int(mCN.group(2))
            TMatrix_contribs = np.empty((rotN,TMatrices.shape[0],2,nelem,2,nelem), dtype=np.complex_)
            for t in targets:
                for i in range(rotN):
                    rotangle = 2*np.pi*i / rotN
                    rot =  qpms.WignerD_yy_fromvector(lMax,np.array([0,0,rotangle]))
                    rotinv = qpms.WignerD_yy_fromvector(lMax,np.array([0,0,-rotangle]))
                    TMatrix_contribs[i] = qpms.apply_matrix_left(rot,qpms.apply_matrix_left(rotinv, TMatrices[:,t], -3),-1)
                TMatrices[:,t] = np.sum(TMatrix_contribs, axis=0) / rotN
        else:
            raise
    elif op[1] == 'tr':
        mCN = reCN.match(op[2]) # Fuck van Rossum for not having assignments inside expressions
        if op[2] == 'σ_z':
            for t in targets:
                TMatrices[:,t] = qpms.apply_ndmatrix_left(zfl,qpms.apply_ndmatrix_left(zfl, TMatrices[:,t], (-4,-3)),(-2,-1))
        elif op[2] == 'σ_y':
            for t in targets:
                TMatrices[:,t] = qpms.apply_ndmatrix_left(yfl,qpms.apply_ndmatrix_left(yfl, TMatrices[:,t], (-4,-3)),(-2,-1))
        elif op[2] == 'σ_x':
            for t in targets:
                TMatrices[:,t] = qpms.apply_ndmatrix_left(xfl,qpms.apply_ndmatrix_left(xfl, TMatrices[:,t], (-4,-3)),(-2,-1))
        elif op[2] == 'C2':
            for t in targets:
                TMatrices[:,t] = qpms.apply_matrix_left(c2rot,qpms.apply_matrix_left(c2rot, TMatrices[:,t], -3),-1)       
        elif mCN:
            rotN = int(mCN.group(2))
            power = int(mCN.group(1)) if mCN.group(1) else 1
            TMatrix_contribs = np.empty((rotN,TMatrices.shape[0],2,nelem,2,nelem), dtype=np.complex_)
            for t in targets:
                rotangle = 2*np.pi*power/rotN
                rot = qpms.WignerD_yy_fromvector(lMax, np.array([0,0,rotangle]))
                rotinv = qpms.WignerD_yy_fromvector(lMax, np.array([0,0,-rotangle]))
                TMatrices[:,t] = qpms.apply_matrix_left(rot, qpms.apply_matrix_left(rotinv, TMatrices[:,t], -3),-1)
        else:
            raise
    elif op[1] == 'copy':
        raise # not implemented
    elif op[1] == 'mult':
        raise # not implemented
    elif op[1] == 'multl':
        incy = np.full((nelem,), False, dtype=bool)
        for incl in op[2][0].split(','):
            l = int(incl)
            incy += (l == ny)
        scaty = np.full((nelem,), False, dtype=bool)
        for scatl in op[2][1].split(','):
            l = int(scatl)
            scaty += (l == ny)
        for t in targets:
            TMatrices[np.ix_(np.arange(TMatrices.shape[0]),np.array([t]),np.array([0,1]),scaty,np.array([0,1]),incy)] *= float(op[2][2])
    else:
        raise #unknown operation; should not happen

TMatrices_interp = interpolate.interp1d(freqs_orig*interpfreqfactor, TMatrices, axis=0, kind='linear',fill_value="extrapolate")

klist_full = qpms.generate_trianglepoints(kdensity, v3d=True, include_origin=True)*3*math.pi/(3*kdensity*hexside)
TMatrices_om = TMatrices_interp(freq)

chunkn = math.ceil(klist_full.shape[0] / chunklen)

if verbose:
    print('Evaluating %d k-points in %d chunks' % (klist_full.shape[0], chunkn), file = sys.stderr)
    sys.stderr.flush()

metadata = np.array({
		'lMax' : lMax,
                'maxlayer' : maxlayer,
                'gaussianSigma' : gaussianSigma,
                'epsilon_b' : epsilon_b,
		'hexside' : hexside,
                'chunkn' : chunkn,
                'TMatrix_file' : TMatrix_file,
                'ops' : ops,
                })

for chunki in range(chunkn):
    svdout = '%s_%dnm_%.4f_c%03d.npz' % (pargs.output_prefix, hexside/1e-9, eVfreq, chunki)

    klist = klist_full[chunki * chunklen : (chunki + 1) * chunklen]

    svdres = qpms.hexlattice_zsym_getSVD(lMax=lMax, TMatrices_om=TMatrices_om, epsilon_b=epsilon_b, hexside=hexside, maxlayer=maxlayer,
            omega=freq, klist=klist, gaussianSigma=gaussianSigma, onlyNmin=False, verbose=verbose)

    #((svUfullTElist, svSfullTElist, svVfullTElist), (svUfullTMlist, svSfullTMlist, svVfullTMlist)) = svdres

    np.savez(svdout, omega = freq, klist = klist,
            metadata=metadata,
                     uTE = svdres[0][0],
                     vTE = svdres[0][2],
                     sTE = svdres[0][1],
                     uTM = svdres[1][0],
                     vTM = svdres[1][2],
                     sTM = svdres[1][1],

    )
    svdres=None

    if scp_dest:
        if svdout:
            subprocess.run(['scp', svdout, scp_dest])

_time_e(btime, verbose)
#print(time.strftime("%H.%M:%S",time.gmtime(time.time()-begtime)))
