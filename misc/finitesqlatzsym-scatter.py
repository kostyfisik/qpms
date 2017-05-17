#!/usr/bin/env python3

import argparse, re, random, string
import subprocess
from scipy.constants import hbar, e as eV, pi, c

unitcell_size = 1 # rectangular lattice
unitcell_indices = tuple(range(unitcell_size))

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
parser.add_argument('--output_prefix', '-p', '-o', action='store', required=True, help='Prefix to the npz output (will be appended frequency, hexside and chunkno)')
parser.add_argument('--nosuffix', action='store_true', help='Do not add dimension metadata to the output filenames')
#sizepar = parser.add_mutually_exclusive_group(required=True)
#parser.add_argument('--hexside', action='store', type=float, required=True, help='Lattice hexagon size length')
parser.add_argument('--dx', action='store', type=float, required=True, help='x-direction lattice constant')
parser.add_argument('--dy', action='store', type=float, required=True, help='y-direction lattice constant')

parser.add_argument('--Nx', '--nx', action='store', type=int, required=True, help='Lattice points in the x-direction')
parser.add_argument('--Ny', '--ny', action='store', type=int, required=True, help='Lattice points in the y-direction')

# In these default settings, the area is 2x2 times larger than first BZ
parser.add_argument('--kxmin', action='store', type=float, default=-1., help='TODO')
parser.add_argument('--kxmax', action='store', type=float, default=1., help='TODO')
parser.add_argument('--kymin', action='store', type=float, default=-1., help='TODO')
parser.add_argument('--kymax', action='store', type=float, default=1., help='TODO')

#parser.add_argument('--kdensity', action='store', type=int, default=33, help='Number of k-points per x-axis segment')
parser.add_argument('--kxdensity', action='store', type=int, default=51, help='k-space resolution in the x-direction')
parser.add_argument('--kydensity', action='store', type=int, default=51, help='k-space resolution in the y-direction')

partgrp = parser.add_mutually_exclusive_group()
partgrp.add_argument('--only_TE', action='store_true', help='Calculate only the projection on the E⟂z modes')
partgrp.add_argument('--only_TM', action='store_true', help='Calculate only the projection on the E∥z modes')
partgrp.add_argument('--serial', action='store_true', help='Calculate the TE and TM parts separately to save memory')

parser.add_argument('--nocentre', action='store_true', help='Place the coordinate origin to the left bottom corner rather that to the centre of the array')

parser.add_argument('--plot_TMatrix', action='store_true', help='Visualise TMatrix on the first page of the output')
#parser.add_argument('--SVD_output', action='store', help='Path to output singular value decomposition result')
parser.add_argument('--maxlayer', action='store', type=int, default=100, help='How far to sum the lattice points to obtain the dispersion')
parser.add_argument('--scp_to', action='store', metavar='N', type=str, help='SCP the output files to a given destination')
parser.add_argument('--background_permittivity', action='store', type=float, default=1., help='Background medium relative permittivity (default 1)')
parser.add_argument('--eVfreq', action='store', required=True, type=float, help='Frequency in eV')

parser.add_argument('--chunklen', action='store', type=int, default=3000, help='Number of k-points per output file (default 3000)')
parser.add_argument('--lMax', action='store', type=int, help='Override lMax from the TMatrix file')
#TODO some more sophisticated x axis definitions
#parser.add_argument('--gaussian', action='store', type=float, metavar='σ', help='Use a gaussian envelope for weighting the interaction matrix contributions (depending on the distance), measured in unit cell lengths (?) FIxME).')
parser.add_argument('--verbose', '-v', action='count', help='Be verbose (about computation times, mostly)')
popgrp=parser.add_argument_group(title='Operations')
popgrp.add_argument('--tr', dest='ops', action=make_action_sharedlist('tr', 'ops'), default=list()) # the default value for dest can be set once
for i in unitcell_indices:
    popgrp.add_argument('--tr%d'%i, dest='ops', action=make_action_sharedlist('tr%d'%i, 'ops'))
popgrp.add_argument('--sym', dest='ops', action=make_action_sharedlist('sym', 'ops'))
for i in unitcell_indices:
    popgrp.add_argument('--sym%d'%i, dest='ops', action=make_action_sharedlist('sym%d'%i, 'ops'))
#popgrp.add_argument('--mult', dest='ops', nargs=3, metavar=('INCSPEC', 'SCATSPEC', 'MULTIPLIER'), action=make_action_sharedlist('mult', 'ops'))
#popgrp.add_argument('--mult0', dest='ops', nargs=3, metavar=('INCSPEC', 'SCATSPEC', 'MULTIPLIER'), action=make_action_sharedlist('mult0', 'ops'))
#popgrp.add_argument('--mult1', dest='ops', nargs=3, metavar=('INCSPEC', 'SCATSPEC', 'MULTIPLIER'), action=make_action_sharedlist('mult1', 'ops'))
popgrp.add_argument('--multl', dest='ops', nargs=3, metavar=('INCL[,INCL,...]', 'SCATL[,SCATL,...]', 'MULTIPLIER'), action=make_action_sharedlist('multl', 'ops'))
for i in unitcell_indices:
    popgrp.add_argument('--multl%d'%i, dest='ops', nargs=3, metavar=('INCL[,INCL,...]', 'SCATL[,SCATL,...]', 'MULTIPLIER'), action=make_action_sharedlist('multl%d'%i, 'ops'))
#popgrp.add_argument('--multl1', dest='ops', nargs=3, metavar=('INCL[,INCL,...]', 'SCATL[,SCATL,...]', 'MULTIPLIER'), action=make_action_sharedlist('multl1', 'ops'))
parser.add_argument('--frequency_multiplier', action='store', type=float, default=1., help='Multiplies the frequencies in the TMatrix file by a given factor.')
# TODO enable more flexible per-sublattice specification
pargs=parser.parse_args()
print(pargs)

maxlayer=pargs.maxlayer
eVfreq = pargs.eVfreq
freq = eVfreq*eV/hbar
verbose=pargs.verbose
dy = pargs.dy
dx = pargs.dx
Ny = pargs.Ny
Nx = pargs.Nx

TMatrix_file = pargs.TMatrix

epsilon_b = pargs.background_permittivity #2.3104
#gaussianSigma = pargs.gaussian if pargs.gaussian else None # hexside * 222 / 7
interpfreqfactor = pargs.frequency_multiplier
scp_dest = pargs.scp_to if pargs.scp_to else None
kxdensity = pargs.kxdensity
kydensity = pargs.kydensity
chunklen = pargs.chunklen

ops = list()
opre = re.compile('(tr|sym|copy|multl|mult)(\d*)')
for oparg in pargs.ops:
    opm = opre.match(oparg[0])
    if opm:
        ops.append(((opm.group(2),) if opm.group(2) else unitcell_indices, opm.group(1), oparg[1]))
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
refind = math.sqrt(epsilon_b)
cdn = c / refind
k_0 = freq * refind / c # = freq / cdn
TMatrices_orig, freqs_orig, freqs_weirdunits_orig, lMaxTM = qpms.loadScuffTMatrices(TMatrix_file)
lMax = lMaxTM
if pargs.lMax:
    lMax = pargs.lMax if pargs.lMax else lMaxTM    
my, ny = qpms.get_mn_y(lMax)
nelem = len(my)
print(TMatrices_orig.shape)
if pargs.lMax: #force commandline specified lMax
    TMatrices_orig = TMatrices_orig[...,0:nelem,:,0:nelem]
print(TMatrices_orig.shape)
TMatrices = np.array(np.broadcast_to(TMatrices_orig[:,nx,:,:,:,:],(len(freqs_orig),unitcell_size,2,nelem,2,nelem)) )
print(TMatrices.shape)
xfl = qpms.xflip_tyty(lMax)
yfl = qpms.yflip_tyty(lMax)
zfl = qpms.zflip_tyty(lMax)
c2rot = qpms.apply_matrix_left(qpms.yflip_yy(3),qpms.xflip_yy(3),-1)

reCN = re.compile('(\d*)C(\d+)')
#TODO C nekonečno

for op in ops:
    if op[0] == 'all':
        #targets = (0,1)
        targets = unitcell_indices
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

print(TMatrices.shape)
TMatrices_interp = interpolate.interp1d(freqs_orig*interpfreqfactor, TMatrices, axis=0, kind='linear',fill_value="extrapolate")


xpositions = np.arange(Nx) * dx
ypositions = np.arange(Ny) * dy
if not pargs.nocentre:
    xpositions -= Nx * dx / 2
    ypositions -= Ny * dy / 2
xpositions, ypositions = np.meshgrid(xpositions, ypositions, indexing='ij', copy=False)
positions=np.stack((xpositions.ravel(),ypositions.ravel()), axis=-1)
N = positions.shape[0]

kx = np.linspace(pargs.kxmin, pargs.kxmax, num=pargs.kxdensity, endpoint=True) * 2*np.pi / dx
ky = np.linspace(pargs.kymin, pargs.kymax, num=pargs.kydensity, endpoint=True) * 2*np.pi / dy
kx, ky = np.meshgrid(kx, ky, indexing='ij', copy=False)
kz = np.sqrt(k_0 - (kx ** 2 + ky ** 2))

klist_full = np.stack((kx,ky,kz), axis=-1).reshape((-1,3))
TMatrices_om = TMatrices_interp(freq)
print(TMatrices_om.shape)

chunkn = math.ceil(klist_full.size / 3 / chunklen)

if verbose:
    print('Evaluating %d k-points in %d chunks' % (klist_full.size / 3, chunkn), file = sys.stderr)
    sys.stderr.flush()

try:
    version = qpms.__version__
except NameError:
    version = None

metadata = np.array({
                'script': os.path.basename(__file__),
                'version': version,
                'type' : 'Plane wave scattering on a finite rectangular lattice',
		'lMax' : lMax,
                'dx' : dx,
                'dy' : dy,
                'Nx' : Nx,
                'Ny' : Ny,
                #'maxlayer' : maxlayer,
                #'gaussianSigma' : gaussianSigma,
                'epsilon_b' : epsilon_b,
		#'hexside' : hexside,
                'chunkn' : chunkn,
                'chunki' : 0,
                'TMatrix_file' : TMatrix_file,
                'ops' : ops,
                'centred' : not pargs.nocentre
                })


scat = qpms.Scattering_2D_zsym(positions, TMatrices_om, k_0, verbose=verbose)

if pargs.only_TE:
    actions = (0,)
elif pargs.only_TM:
    actions = (1,)
elif pargs.serial:
    actions = (0,1)
else:
    actions = (None,)

xu = np.array((1,0,0))
yu = np.array((0,1,0))
zu = np.array((0,0,1))
TEč, TMč = qpms.symz_indexarrays(lMax)

klist_full_2D = klist_full[...,:2]
klist_full_dir = klist_full/np.linalg.norm(klist_full, axis=-1, keepdims=True)
for action in actions:
    if action is None:
        scat.prepare(verbose=verbose)
        actionstring = ''
    else:
        scat.prepare_partial(action, verbose=verbose)
        actionstring = '.TM' if action else '.TE'
    for chunki in range(chunkn):
        if pargs.nosuffix:
            outfile = pargs.output_prefix + actionstring + (
                    ('.%03d' % chunki) if chunkn > 1 else '')
        else:
            outfile = '%s_%dx%d_%.0fnmx%.0fnm_%.4f%s%s.npz' % (
                    pargs.output_prefix, Nx, Ny, dx/1e-9, dy/1e-9,
                    eVfreq, actionstring,
                    (".%03d" % cunki) if chunkn > 1 else '')

        klist = klist_full[chunki * chunklen : (chunki + 1) * chunklen]
        klist2d = klist_full_2D[chunki * chunklen : (chunki + 1) * chunklen]
        klistdir = klist_full_dir[chunki * chunklen : (chunki + 1) * chunklen]
        
        '''
        The following loop is a fuckup that has its roots in the fact that
        the function qpms.get_π̃τ̃_y1 in qpms_p.py is not vectorized
        (and consequently, neither is plane_pq_y.)

        And Scattering_2D_zsym.scatter_partial is not vectorized, either.
        '''
        if action == 0 or action is None:
            xresult = np.full((klist.shape[0], N, nelem), np.nan, dtype=complex)
            yresult = np.full((klist.shape[0], N, nelem), np.nan, dtype=complex)
        if action == 1 or action is None:
            zresult = np.full((klist.shape[0], N, nelem), np.nan, dtype=complex)
        for i in range(klist.shape[0]):
            if math.isnan(klist[i,2]):
                continue
            kdir = klistdir[i]
            phases = np.exp(np.sum(klist2d[i] * positions, axis=-1))
            if action == 0 or action is None:
                pq = np.array(qpms.plane_pq_y(lMax, kdir, xu)).ravel()[TEč] * phases[:, nx] 
                xresult[i] = scat.scatter_partial(0, pq)
                pq = np.array(qpms.plane_pq_y(lMax, kdir, yu)).ravel()[TEč] * phases[:, nx] 
                yresult[i] = scat.scatter_partial(0, pq)
            if action == 1 or action is None:
                pq = np.array(qpms.plane_pq_y(lMax, kdir, xu)).ravel()[TMč] * phases[:, nx] 
                zresult[i] = scat.scatter_partial(1, pq)

        metadata[()]['chunki'] = chunki
        if action is None:
            np.savez(outfile, omega = freq, klist = klist,
                metadata=metadata,
                ab_x=xresult,
                ab_y=yresult,
                ab_z=zresult
                )
        elif action == 0:
            np.savez(outfile, omega = freq, klist = klist,
                metadata=metadata,
                ab_x=xresult,
                ab_y=yresult,
                )
        elif action == 1:
            np.savez(outfile, omega = freq, klist = klist,
                metadata=metadata,
                ab_z=zresult
                )
        else:
            raise

        if scp_dest:
            if outfile:
                subprocess.run(['scp', outfile, scp_dest])
    scat.forget_matrices() # free memory in case --serial was used

_time_e(btime, verbose)
