#!/usr/bin/env python3

import argparse, re, random, string
from scipy.constants import hbar, e as eV, pi, c

parser = argparse.ArgumentParser()
#TODO? použít type=argparse.FileType('r') ?
parser.add_argument('--output', '-o', action='store', help='Path to output PDF')
parser.add_argument('--nSV', action='store', metavar='N', type=int, default=1, help='Draw N minimum singular values')
parser.add_argument('--bitmap', action='store_true', help='Create an interpolated bitmap rather than a scatter plot.')
#parser.add_argument('--eVfreq', action='store', required=True, type=float, help='Frequency in eV')
parser.add_argument('inputfile',  nargs='+', help='Npz file(s) generated by dispersion_chunks.py or other script')
pargs=parser.parse_args()
print(pargs)

#freq = eVfreq*eV/hbar

pdfout = pargs.output if pargs.output else '%s.pdf' % pargs.inputfile[-1]
print(pdfout)

svn = pargs.nSV

# -----------------finished basic CLI parsing (except for op arguments) ------------------
import time
begtime=time.time()

import qpms
from matplotlib.path import Path
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np
import os, sys, warnings, math
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy import interpolate

# We do not want to import whole qpms, so copy and modify this only fun needed
def nelem2lMax(nelem):
    lMax = round(math.sqrt(1+nelem) - 1)
    if ((lMax < 1) or ((lMax + 2) * lMax != nelem)):
        raise
    else:
        return lMax



nx = None
s3 = math.sqrt(3)


# read data from files
lMax = None
epsilon_b = None
hexside = None
karrlist = list()
svTElist = list()
svTMlist = list()
omegalist = list()
records = 0
for filename in pargs.inputfile:
    npz = np.load(filename)
    lMaxRead = npz['metadata'][()]['lMax'] if 'lMax' in npz['metadata'][()] else nelem2lMax(npz['sTE'].shape[1] / 2)
    if lMax is None: lMax = lMaxRead
    elif lMax != lMaxRead: raise
    if epsilon_b is None: epsilon_b = npz['metadata'][()]['epsilon_b']
    elif epsilon_b != npz['metadata'][()]['epsilon_b'] : raise
    if hexside is None: hexside = npz['metadata'][()]['hexside']
    elif hexside != npz['metadata'][()]['hexside'] : raise
    omegalist.append(npz['omega'][()])
    karrlist.append(np.array(npz['klist']))
    svTElist.append(np.array(npz['sTE'][:,-svn:]))
    svTMlist.append(np.array(npz['sTM'][:,-svn:]))
    records += 1
    npz.close()

# sort by frequencies
omegas = set(omegalist)
print(omegas)
k = dict()
svTE = dict()
svTM = dict()
for omega in omegas:
    k[omega] = list()
    svTE[omega] = list()
    svTM[omega] = list()
for i in range(records):
    omega = omegalist[i]
    k[omega].append(karrlist[i])
    svTE[omega].append(svTElist[i])
    svTM[omega].append(svTMlist[i])
# concatenate arrays for each frequency
for omega in omegas:
    k[omega] = np.concatenate(k[omega])
    svTE[omega] = np.concatenate(svTE[omega])
    svTM[omega] = np.concatenate(svTM[omega])

# ... that was for the slices. TODO fill also the righternmost plot with the calculated (which?) modes.

pdf = PdfPages(pdfout)

# In[3]:

cdn = c/ math.sqrt(epsilon_b)
#my, ny = qpms.get_mn_y(lMax)
#nelem = len(my)
nelem = lMax * (lMax + 2)


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


for omega in sorted(omegas):
    klist = k[omega]
    if pargs.bitmap:
            minx = np.amin(klist[:,0])
            maxx = np.amax(klist[:,0])
            miny = np.amin(klist[:,1])
            maxy = np.amax(klist[:,1])
            l = klist.shape[0]
            meshstep = math.sqrt((maxy - miny) * (maxx - minx) / l) / 9
            x = np.linspace(minx, maxx, (maxx-minx) / meshstep)
            y = np.linspace(miny, maxy, (maxy-miny) / meshstep)
            fullshape = np.broadcast(x[:,nx],y[nx,:]).shape
            flatx = np.broadcast_to(x[:,nx], fullshape).flatten
            flaty = np.broadcast_to(y[nx,:], fullshape).flatten
    minsvTElist = svTE[omega]
    minsvTMlist = svTM[omega]
    for minN in reversed(range(svn)):
        f, axes = plt.subplots(1,3, figsize=(20,4.8))
        ax = axes[0]
        if pargs.bitmap:
            interpolator = interpolate.interp2d(klist[:,0], klist[:,1], np.abs(minsvTElist[:,minN]))
            z = interpolator(flatx, flaty)
            z.reshape(fullshape)
            sc = ax.pcolormesh(x[:,nx],y[nx,:],z)
        else:
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
        ax.title.set_text('E in-plane ("TE"), %d. lowest SV' % minN)
        f.colorbar(sc,ax=ax)
        
        
        ax = axes[1]
        if pargs.bitmap:
            interpolator = interpolate.interp2d(klist[:,0], klist[:,1], np.abs(minsvTMlist[:,minN]))
            meshstep = math.sqrt((maxy - miny) * (maxx - minx) / l) / 9
            z = interpolator(flatx, flaty)
            z.reshape(fullshape)
            sc = ax.pcolormesh(x[:,nx],y[nx,:],z)
        else:
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
        ax.title.set_text('E perpendicular ("TM"), %d. lowest SV' % minN)
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

