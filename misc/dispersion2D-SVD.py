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
parser.add_argument('--output_prefix', action='store', required=True, help='Prefix to the pdf and/or npz output (will be appended frequency and hexside)')
#sizepar = parser.add_mutually_exclusive_group(required=True)
parser.add_argument('--hexside', action='store', type=float, required=True, help='Lattice hexagon size length')
parser.add_argument('--output', action='store', help='Path to output PDF')
parser.add_argument('--store_SVD', action='store_false', help='If specified without --SVD_output, it will save the data in a file named as the PDF output, but with .npz extension instead')
parser.add_argument('--plot_TMatrix', action='store_true', help='Visualise TMatrix on the first page of the output')
#parser.add_argument('--SVD_output', action='store', help='Path to output singular value decomposition result')
parser.add_argument('--nSV', action='store', metavar='N', type=int, default=1, help='Store and draw N minimun singular values')
parser.add_argument('--maxlayer', action='store', type=int, default=100, help='How far to sum the lattice points to obtain the dispersion')
parser.add_argument('--scp_to', action='store', metavar='N', type=str, help='SCP the output files to a given destination')
parser.add_argument('--background_permittivity', action='store', type=float, default=1., help='Background medium relative permittivity (default 1)')
parser.add_argument('--eVfreq', action='store', required=True, type=float, help='Frequency in eV')
parser.add_argument('--kdensity', action='store', type=int, default=33, help='Number of k-points per x-axis segment')
parser.add_argument('--lMax', action='store', type=int, help='Override lMax from the TMatrix file')
#TODO some more sophisticated x axis definitions
parser.add_argument('--gaussian', action='store', type=float, metavar='σ', help='Use a gaussian envelope for weighting the interaction matrix contributions (depending on the distance), measured in unit cell lengths (?) FIxME).')
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

TMatrix_file = pargs.TMatrix
pdfout = pargs.output if pargs.output else (
	'%s_%dnm_%.4f.pdf' % (pargs.output_prefix,hexside/1e-9,eVfreq) if pargs.output_prefix else
	(''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(10)) + '.pdf'))
print(pdfout)
if(pargs.store_SVD):
    if re.search('.pdf$', pdfout):
        svdout = re.sub('.pdf$', r'.npz', pdfout)
    else:
        svdout = pdfout + '.npz'
else:
     svdout = None


epsilon_b = pargs.background_permittivity #2.3104
gaussianSigma = pargs.gaussian if pargs.gaussian else None # hexside * 222 / 7
interpfreqfactor = pargs.frequency_multiplier
scp_dest = pargs.scp_to if pargs.scp_to else None
kdensity = pargs.kdensity
svn = pargs.nSV

# TODO multiplier operation definitions and parsing
#factor13inc = 10
#factor13scat=10

ops = list()
opre = re.compile('(tr|sym|copy|multl|mult)(\d*)')
for oparg in pargs.ops:
    opm = opre.match(oparg[0])
    if opm:
        ops.append(((opm.group(2),) if opm.group(2) else (0,1), opm.group(1), oparg[1]))
    else:
        raise # should not happen
print(ops)

#ops = (
#    # co, typ operace (symetrizace / transformace / kopie), specifikace (operace nebo zdroj), 
#    # co: 0, 1, (0,1), (0,), (1,), #NI: 'all'
#    # typ operace: sym, tr, copy
#    # specifikace:
#    # sym, tr: 'σ_z', 'σ_y', 'C2'; sym: 'C3', 
#    # copy: 0, 1 (zdroj)
#    ((0,1), 'sym', 'σ_z'),
#    #((0,1), 'sym', 'σ_x'),
#    #((0,1), 'sym', 'σ_y'),
#    ((0,1), 'sym', 'C3'),
#    ((1), 'tr', 'C2'),
#
#)

# -----------------finished basic CLI parsing (except for op arguments) ------------------
import time
begtime=time.time()

from matplotlib.path import Path
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import qpms
import numpy as np
import os, sys, warnings, math
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy import interpolate
nx = None
s3 = math.sqrt(3)



pdf = PdfPages(pdfout)

# In[3]:

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

ž = np.arange(2*nelem)
tž = ž // nelem
mž = my[ž%nelem]
nž = ny[ž%nelem]
TEž = ž[(mž+nž+tž) % 2 == 0]
TMž = ž[(mž+nž+tž) % 2 == 1]

č = np.arange(2*2*nelem)
žč = č % (2* nelem)
tč = tž[žč]
mč = mž[žč]
nč = nž[žč]
TEč = č[(mč+nč+tč) % 2 == 0]
TMč = č[(mč+nč+tč) % 2 == 1]

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



# In[4]:
if(pargs.plot_TMatrix):    
    om = np.linspace(np.min(freqs_orig), np.max(freqs_orig),100)
    TMatrix0ip = np.reshape(TMatrices_interp(om)[:,0], (len(om), 2*nelem*2*nelem))
    f, axa = plt.subplots(2, 2, figsize=(15,15))
    
    #print(TMatrices.shape)
    #plt.plot(om, TMatrices[:,0,0,0,0].imag,'r',om, TMatrices[:,0,0,0,0].real,'r--',om, TMatrices[:,0,2,0,2].imag,'b',om, TMatrices[:,0,2,0,2].real,'b--'))
            
    ax = axa[0,0]
    ax2 = ax.twiny()
    ax2.set_xlim([ax.get_xlim()[0]/eV*hbar,ax.get_xlim()[1]/eV*hbar])
    ax.plot(
             om, TMatrix0ip[:,:].imag,'-',om, TMatrix0ip[:,:].real,'--',
    )
    ax = axa[0,1]
    ax2 = ax.twiny()
    ax2.set_xlim([ax.get_xlim()[0]/eV*hbar,ax.get_xlim()[1]/eV*hbar])
    ax.plot(
             om, abs(TMatrix0ip[:,:]),'-'
    )
    ax.set_yscale('log')
    
    ax = axa[1,1]
    ax2 = ax.twiny()
    ax2.set_xlim([ax.get_xlim()[0]/eV*hbar,ax.get_xlim()[1]/eV*hbar])
    ax.plot(
             om, np.unwrap(np.angle(TMatrix0ip[:,:]),axis=0),'-'
    )

    ax = axa[1,0]
    ax.text(0.5,0.5,str(pargs).replace(',',',\n'),horizontalalignment='center',verticalalignment='center',transform=ax.transAxes)
    pdf.savefig(f)


# In[ ]:

'''
#kdensity = 66 #defined from cl arguments
bz_0 = np.array((0,0,0.,))
bz_K1 = np.array((1.,0,0))*4*np.pi/3/hexside/s3
bz_K2 = np.array((1./2.,s3/2,0))*4*np.pi/3/hexside/s3
bz_M = np.array((3./4, s3/4,0))*4*np.pi/3/hexside/s3
k0Mlist = bz_0 + (bz_M-bz_0) * np.linspace(0,1,kdensity)[:,nx]
kMK1list = bz_M + (bz_K1-bz_M) * np.linspace(0,1,kdensity)[:,nx]
kK10list = bz_K1 + (bz_0-bz_K1) * np.linspace(0,1,kdensity)[:,nx]
k0K2list = bz_0 + (bz_K2-bz_0) * np.linspace(0,1,kdensity)[:,nx]
kK2Mlist = bz_K2 + (bz_M-bz_K2) * np.linspace(0,1,kdensity)[:,nx]
B1 = 2* bz_K1 - bz_K2
B2 = 2* bz_K2 - bz_K1
klist = np.concatenate((k0Mlist,kMK1list,kK10list,k0K2list,kK2Mlist), axis=0)
kxmaplist = np.concatenate((np.array([0]),np.cumsum(np.linalg.norm(np.diff(klist, axis=0), axis=-1))))
'''
klist = qpms.generate_trianglepoints(kdensity, v3d=True, include_origin=True)*3*math.pi/(3*kdensity*hexside)

# In[ ]:

n2id = np.identity(2*nelem)
n2id.shape = (2,nelem,2,nelem)
nan = float('nan')
filecount = 0
for one in (1,):
    omega = freq # from args; k_0 * c / math.sqrt(epsilon_b)
    k_0 = omega * math.sqrt(epsilon_b) / c
    tdic = qpms.hexlattice_get_AB(lMax,k_0*hexside,maxlayer)
    #print(filecount, omega/eV*hbar)
    #sys.stdout.flush()
    a_self = tdic['a_self'][:,:nelem,:nelem]
    b_self = tdic['b_self'][:,:nelem,:nelem]
    a_u2d = tdic['a_u2d'][:,:nelem,:nelem]
    b_u2d = tdic['b_u2d'][:,:nelem,:nelem]
    a_d2u = tdic['a_d2u'][:,:nelem,:nelem]
    b_d2u = tdic['b_d2u'][:,:nelem,:nelem]
    unitcell_translations = tdic['self_tr']*hexside*s3
    u2d_translations = tdic['u2d_tr']*hexside*s3
    d2u_translations = tdic['d2u_tr']*hexside*s3
    if gaussianSigma:
        unitcell_envelope = np.exp(-np.sum(tdic['self_tr']**2,axis=-1)/(2*gaussianSigma**2))
        u2d_envelope = np.exp(-np.sum(tdic['u2d_tr']**2,axis=-1)/(2*gaussianSigma**2))
        d2u_envelope = np.exp(-np.sum(tdic['d2u_tr']**2,axis=-1)/(2*gaussianSigma**2))
    
    
    TMatrices_om = TMatrices_interp(omega)
    if svdout:
        svUfullTElist = np.full((klist.shape[0], 2*nelem, 2*nelem), np.nan, dtype=complex)
        svVfullTElist = np.full((klist.shape[0], 2*nelem, 2*nelem), np.nan, dtype=complex)
        svSfullTElist = np.full((klist.shape[0], 2*nelem), np.nan, dtype=complex)
        svUfullTMlist = np.full((klist.shape[0], 2*nelem, 2*nelem), np.nan, dtype=complex)
        svVfullTMlist = np.full((klist.shape[0], 2*nelem, 2*nelem), np.nan, dtype=complex)
        svSfullTMlist = np.full((klist.shape[0], 2*nelem), np.nan, dtype=complex)

     
    minsvTElist = np.full((klist.shape[0], svn),np.nan)
    minsvTMlist = np.full((klist.shape[0], svn),np.nan)
    leftmatrixlist = np.full((klist.shape[0],2,2,nelem,2,2,nelem),np.nan,dtype=complex)
    isNaNlist = np.zeros((klist.shape[0]), dtype=bool)
    
    # sem nějaká rozumná smyčka
    for ki in range(klist.shape[0]):
        k = klist[ki]
        if (k_0*k_0 - k[0]*k[0] - k[1]*k[1] < 0):
            isNaNlist[ki] = True
            continue

        phases_self = np.exp(1j*np.tensordot(k,unitcell_translations,axes=(0,-1)))
        phases_u2d = np.exp(1j*np.tensordot(k,u2d_translations,axes=(0,-1)))
        phases_d2u = np.exp(1j*np.tensordot(k,d2u_translations,axes=(0,-1)))
        if gaussianSigma:
            phases_self *= unitcell_envelope
            phases_u2d *= u2d_envelope
            phases_d2u *= d2u_envelope
        leftmatrix = np.zeros((2,2,nelem, 2,2,nelem), dtype=complex)
        #       0: u,E<--u,E
        #       1: d,M<--d,M
        leftmatrix[0,0,:,0,0,:] = np.tensordot(a_self,phases_self, axes=(0,-1)) # u2u, E2E
        leftmatrix[1,0,:,1,0,:] = leftmatrix[0,0,:,0,0,:] # d2d, E2E
        leftmatrix[0,1,:,0,1,:] = leftmatrix[0,0,:,0,0,:] # u2u, M2M
        leftmatrix[1,1,:,1,1,:] = leftmatrix[0,0,:,0,0,:] # d2d, M2M
        leftmatrix[0,0,:,0,1,:] = np.tensordot(b_self,phases_self, axes=(0,-1)) # u2u, M2E
        leftmatrix[0,1,:,0,0,:] = leftmatrix[0,0,:,0,1,:] # u2u, E2M
        leftmatrix[1,1,:,1,0,:] = leftmatrix[0,0,:,0,1,:] # d2d, E2M
        leftmatrix[1,0,:,1,1,:] = leftmatrix[0,0,:,0,1,:] # d2d, M2E
        leftmatrix[0,0,:,1,0,:] = np.tensordot(a_d2u, phases_d2u,axes=(0,-1)) #d2u,E2E
        leftmatrix[0,1,:,1,1,:] = leftmatrix[0,0,:,1,0,:] #d2u, M2M
        leftmatrix[1,0,:,0,0,:] = np.tensordot(a_u2d, phases_u2d,axes=(0,-1)) #u2d,E2E
        leftmatrix[1,1,:,0,1,:] = leftmatrix[1,0,:,0,0,:] #u2d, M2M
        leftmatrix[0,0,:,1,1,:] = np.tensordot(b_d2u, phases_d2u,axes=(0,-1)) #d2u,M2E
        leftmatrix[0,1,:,1,0,:] = leftmatrix[0,0,:,1,1,:] #d2u, E2M
        leftmatrix[1,0,:,0,1,:] = np.tensordot(b_u2d, phases_u2d,axes=(0,-1)) #u2d,M2E
        leftmatrix[1,1,:,0,0,:] = leftmatrix[1,0,:,0,1,:] #u2d, E2M
        #leftmatrix is now the translation matrix T
        for j in range(2):
            leftmatrix[j] = -np.tensordot(TMatrices_om[j], leftmatrix[j], axes=([-2,-1],[0,1]))
            # at this point, jth row of leftmatrix is that of -MT
            leftmatrix[j,:,:,j,:,:] += n2id
        #now we are done, 1-MT

        leftmatrixlist[ki] = leftmatrix


    nnlist = np.logical_not(isNaNlist)
    leftmatrixlist_s = np.reshape(leftmatrixlist,(klist.shape[0], 2*2*nelem,2*2*nelem))[nnlist]
    leftmatrixlist_TE = leftmatrixlist_s[np.ix_(np.arange(leftmatrixlist_s.shape[0]),TEč,TEč)]
    leftmatrixlist_TM = leftmatrixlist_s[np.ix_(np.arange(leftmatrixlist_s.shape[0]),TMč,TMč)]
    #svarr = np.linalg.svd(leftmatrixlist_TE, compute_uv=False)
    #argsortlist = np.argsort(svarr, axis=-1)[...,:svn]
    #minsvTElist[nnlist] = svarr[...,argsortlist]
    #minsvTElist[nnlist] = np.amin(np.linalg.svd(leftmatrixlist_TE, compute_uv=False), axis=-1)
    if svdout:
        svUfullTElist[nnlist], svSfullTElist[nnlist], svVfullTElist[nnlist] = np.linalg.svd(leftmatrixlist_TE, compute_uv=True)
        svUfullTMlist[nnlist], svSfullTMlist[nnlist], svVfullTMlist[nnlist] = np.linalg.svd(leftmatrixlist_TM, compute_uv=True)
    minsvTElist[nnlist] = np.linalg.svd(leftmatrixlist_TE, compute_uv=False)[...,-svn:]
    #svarr = np.linalg.svd(leftmatrixlist_TM, compute_uv=False)
    #argsortlist = np.argsort(svarr, axis=-1)[...,:svn]
    #minsvTMlist[nnlist] = svarr[...,argsortlist]
    #minsvTMlist[nnlist] = np.amin(np.linalg.svd(leftmatrixlist_TM, compute_uv=False), axis=-1)
    minsvTMlist[nnlist] = np.linalg.svd(leftmatrixlist_TM, compute_uv=False)[...,-svn:]

                    
'''
omlist = np.broadcast_to(omegalist[:,nx], minsvTElistarr[...,0].shape)
kxmlarr = np.broadcast_to(kxmaplist[nx,:], minsvTElistarr[...,0].shape)
klist = np.concatenate((k0Mlist,kMK1list,kK10list,k0K2list,kK2Mlist), axis=0)
'''

''' The new pretty diffracted order drawing '''
maxlayer_reciprocal=4 
cdn = c/ math.sqrt(epsilon_b)
bz_0 = np.array((0,0,))
bz_K1 = np.array((1.,0))*4*np.pi/3/hexside/s3
bz_K2 = np.array((1./2.,s3/2))*4*np.pi/3/hexside/s3
bz_M = np.array((3./4, s3/4))*4*np.pi/3/hexside/s3

# reciprocal lattice basis
B1 = 2* bz_K1 - bz_K2
B2 = 2* bz_K2 - bz_K1

if svdout:
    np.savez(svdout, omega = freq, klist = klist, bzpoints = np.array([bz_0, bz_K1, bz_K2, bz_M, B1, B2]), 
                     uTE = svUfullTElist,
                     vTE = svVfullTElist,
                     sTE = svSfullTElist,
                     uTM = svUfullTMlist,
                     vTM = svVfullTMlist,
                     sTM = svSfullTMlist,
    )
 
k2density = 100
k0Mlist = bz_0 + (bz_M-bz_0) * np.linspace(0,1,k2density)[:,nx]
kMK1list = bz_M + (bz_K1-bz_M) * np.linspace(0,1,k2density)[:,nx]
kK10list = bz_K1 + (bz_0-bz_K1) * np.linspace(0,1,k2density)[:,nx]
k0K2list = bz_0 + (bz_K2-bz_0) * np.linspace(0,1,k2density)[:,nx]
kK2Mlist = bz_K2 + (bz_M-bz_K2) * np.linspace(0,1,k2density)[:,nx]
k2list = np.concatenate((k0Mlist,kMK1list,kK10list,k0K2list,kK2Mlist), axis=0)
kxmaplist = np.concatenate((np.array([0]),np.cumsum(np.linalg.norm(np.diff(k2list, axis=0), axis=-1))))

centers2=qpms.generate_trianglepoints(maxlayer_reciprocal, v3d = False, include_origin= True)*4*np.pi/3/hexside
rot90 = np.array([[0,-1],[1,0]])
centers2=np.dot(centers2,rot90)

import matplotlib.pyplot as plt
import matplotlib
from matplotlib.path import Path
import matplotlib.patches as patches
cmap = matplotlib.cm.prism
colormax = np.amax(np.linalg.norm(centers2,axis=0))


# In[ ]:
for minN in reversed(range(svn)):
    f, axes = plt.subplots(1,3, figsize=(20,4.8))
    ax = axes[0]
    sc = ax.scatter(klist[:,0], klist[:,1], c = np.clip(np.abs(minsvTElist[:,minN]),0,1), lw=0)
    for center in centers2:
        circle=plt.Circle((center[0],center[1]),omega/cdn, facecolor='none', edgecolor=cmap(np.linalg.norm(center)/colormax),lw=0.5)
        ax.add_artist(circle)
    verts = [(math.cos(math.pi*i/3)*4*np.pi/3/hexside/s3,math.sin(math.pi*i/3)*4*np.pi/3/hexside/s3) for i in range(6 +1)]
    codes = [Path.MOVETO,Path.LINETO,Path.LINETO,Path.LINETO,Path.LINETO,Path.LINETO,Path.CLOSEPOLY,]
    path = Path(verts, codes)
    patch = patches.PathPatch(path, facecolor='none', edgecolor='black',  lw=1)
    ax.add_patch(patch)
    ax.set_xticks([]) 
    ax.set_yticks([]) 
    ax.title.set_text('E in-plane ("TE")')
    f.colorbar(sc,ax=ax)
    
    
    ax = axes[1]
    sc = ax.scatter(klist[:,0], klist[:,1], c = np.clip(np.abs(minsvTMlist[:,minN]),0,1), lw=0)
    for center in centers2:
        circle=plt.Circle((center[0],center[1]),omega/cdn, facecolor='none', edgecolor=cmap(np.linalg.norm(center)/colormax),lw=0.5)
        ax.add_artist(circle)
    verts = [(math.cos(math.pi*i/3)*4*np.pi/3/hexside/s3,math.sin(math.pi*i/3)*4*np.pi/3/hexside/s3) for i in range(6 +1)]
    codes = [Path.MOVETO,Path.LINETO,Path.LINETO,Path.LINETO,Path.LINETO,Path.LINETO,Path.CLOSEPOLY,]
    path = Path(verts, codes)
    patch = patches.PathPatch(path, facecolor='none', edgecolor='black',  lw=1)
    ax.add_patch(patch)
    ax.set_xticks([]) 
    ax.set_yticks([])
    ax.title.set_text('E perpendicular ("TM")')
    f.colorbar(sc,ax=ax)

    ax = axes[2]
    for center in centers2:
        ax.plot(kxmaplist, np.linalg.norm(k2list-center,axis=-1)*cdn, '-', color=cmap(np.linalg.norm(center)/colormax))

    #ax.set_xlim([np.min(kxmlarr),np.max(kxmlarr)])
    #ax.set_ylim([np.min(omegalist),np.max(omegalist)])
    xticklist = [0, kxmaplist[len(k0Mlist)-1], kxmaplist[len(k0Mlist)+len(kMK1list)-1], kxmaplist[len(k0Mlist)+len(kMK1list)+len(kK10list)-1], kxmaplist[len(k0Mlist)+len(kMK1list)+len(kK10list)+len(k0K2list)-1], kxmaplist[len(k0Mlist)+len(kMK1list)+len(kK10list)+len(k0K2list)+len(kK2Mlist)-1]]
    ax.set_xticks(xticklist)
    for xt in xticklist:
        ax.axvline(xt, ls='dotted', lw=0.3,c='k')
    ax.set_xticklabels(['Γ', 'M', 'K', 'Γ', 'K\'','M'])
    ax.axhline(omega, c='black')
    ax.set_ylim([0,5e15])
    ax2 = ax.twinx()
    ax2.set_ylim([ax.get_ylim()[0]/eV*hbar,ax.get_ylim()[1]/eV*hbar])
   
    pdf.savefig(f)
    
pdf.close()

if scp_dest:
    subprocess.run(['scp', pdfout, scp_dest])
    if svdout:
        subprocess.run(['scp', svdout, scp_dest])

print(time.strftime("%H.%M:%S",time.gmtime(time.time()-begtime)))
