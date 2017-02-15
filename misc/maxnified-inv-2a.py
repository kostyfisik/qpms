
# coding: utf-8

# In[1]:

translations_dir = '/l/necadam1/translations-precalc/diracpoints-newdata/222'
TMatrix_file ='/m/home/home4/46/necadam1/unix/tmatrix-experiments/twisted_triangle/silver/twisted_triangle.TMatrix.nonan'

pdfout = '/m/home/home4/46/necadam1/unix/tmp/pdf_out/inv-2-mag10-10.pdf'

hexside = 375e-9
epsilon_b = 2.3104
gaussianSigma = None # hexside * 222 / 7

factor13inc = 10
factor13scat=10

ops = (
    # co, typ operace (symetrizace / transformace / kopie), specifikace (operace nebo zdroj), 
    # co: 0, 1, (0,1), (0,), (1,), #NI: 'all'
    # typ operace: sym, tr, copy
    # specifikace:
    # sym, tr: 'σ_z', 'σ_y', 'C2'; sym: 'C3', 
    # copy: 0, 1 (zdroj)
    ((0,1), 'sym', 'σ_z'),
    #((0,1), 'sym', 'σ_x'),
    #((0,1), 'sym', 'σ_y'),
    ((0,1), 'sym', 'C3'),
    ((1), 'tr', 'C2'),

)

interpfreqfactor = 0.5


import qpms
import numpy as np
import os, sys
import warnings
import math
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.constants import hbar, e as eV, pi, c
from scipy import interpolate
nx = None
s3 = math.sqrt(3)

pdf = PdfPages(pdfout)


# In[2]:

#TODO později
#import argparse
#parser = argparse.ArgumentParser()
#parser.add_argument('--sym', 'mz', 'my', 'mx', 'C3', 'C2' type=str, help='symmetrize both particles')
#args = parser.parse_args()


# In[3]:

# specifikace T-matice zde
cdn = c/ math.sqrt(epsilon_b)
TMatrices_orig, freqs_orig, freqs_weirdunits_orig, lMax = qpms.loadScuffTMatrices(TMatrix_file)
my, ny = qpms.get_mn_y(lMax)
nelem = len(my)

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

TMatrices[:,:,:,:,:,ny==3] *= factor13inc
TMatrices[:,:,:,ny==3,:,:] *= factor13scat
xfl = qpms.xflip_tyty(lMax)
yfl = qpms.yflip_tyty(lMax)
zfl = qpms.zflip_tyty(lMax)
c2rot = qpms.apply_matrix_left(qpms.yflip_yy(3),qpms.xflip_yy(3),-1)

for op in ops:
    if op[0] == 'all':
        targets = (0,1)
    elif isinstance(op[0],int):
        targets = (op[0],)
    else:
        targets = op[0]
        
    if op[1] == 'sym':
        if op[2] == 'σ_z':
            for t in targets:
                TMatrices[:,t] = (TMatrices[:,t] + qpms.apply_ndmatrix_left(zfl,qpms.apply_ndmatrix_left(zfl, TMatrices[:,t], (-4,-3)),(-2,-1)))/2
        elif op[2] == 'σ_y':
            for t in targets:
                TMatrices[:,t] = (TMatrices[:,t] + qpms.apply_ndmatrix_left(yfl,qpms.apply_ndmatrix_left(yfl, TMatrices[:,t], (-4,-3)),(-2,-1)))/2
        elif op[2] == 'σ_x':
            for t in targets:
                TMatrices[:,t] = (TMatrices[:,t] + qpms.apply_ndmatrix_left(xfl,qpms.apply_ndmatrix_left(xfl, TMatrices[:,t], (-4,-3)),(-2,-1)))/2
        elif op[2] == 'C3': # FIXME fuj fuj fuj, použij regex!!!
            rotN = 3
            TMatrix_contribs = np.empty((rotN,TMatrices.shape[0],2,nelem,2,nelem), dtype=np.complex_)
            for t in targets:
                for i in range(rotN):
                    rotangle = 2*np.pi*i / rotN
                    rot =  qpms.WignerD_yy_fromvector(lMax,np.array([0,0,rotangle]))
                    rotinv = qpms.WignerD_yy_fromvector(lMax,np.array([0,0,-rotangle]))
                    TMatrix_contribs[i] = qpms.apply_matrix_left(rot,qpms.apply_matrix_left(rotinv, TMatrices[:,t], -3),-1)
                TMatrices[:,t] = np.sum(TMatrix_contribs, axis=0) / rotN
        elif op[2] == 'C2':
            for t in targets:
                TMatrices[:,t] = (TMatrices[:,t] + qpms.apply_matrix_left(c2rot,qpms.apply_matrix_left(c2rot, TMatrices[:,t], -3),-1))/2                
        else:
            raise
    elif op[1] == 'tr':
        if op[2] == 'σ_z':
            for t in targets:
                TMatrices[:,t] = qpms.apply_ndmatrix_left(zfl,qpms.apply_ndmatrix_left(zfl, TMatrices[:,t], (-4,-3)),(-2,-1))
        elif op[2] == 'σ_y':
            for t in targets:
                TMatrices[:,t] = qpms.apply_ndmatrix_left(yfl,qpms.apply_ndmatrix_left(yfl, TMatrices[:,t], (-4,-3)),(-2,-1))
        elif op[2] == 'σ_x':
            for t in targets:
                TMatrices[:,t] = qpms.apply_ndmatrix_left(xfl,qpms.apply_ndmatrix_left(xfl, TMatrices[:,t], (-4,-3)),(-2,-1))
        elif op[2] == 'C3': # TODO use regex and generalize
            rotN = 3
            TMatrix_contribs = np.empty((rotN,TMatrices.shape[0],2,nelem,2,nelem), dtype=np.complex_)
            for t in targets:
                for i in range(rotN):
                    rotangle = 2*np.pi*i / rotN
                    rot =  qpms.WignerD_yy_fromvector(lMax,np.array([0,0,rotangle]))
                    rotinv = qpms.WignerD_yy_fromvector(lMax,np.array([0,0,-rotangle]))
                    TMatrix_contribs[i] = qpms.apply_matrix_left(rot,qpms.apply_matrix_left(rotinv, TMatrices[:,t], -3),-1)
        elif op[2] == 'C2':
            for t in targets:
                TMatrices[:,t] = qpms.apply_matrix_left(c2rot,qpms.apply_matrix_left(c2rot, TMatrices[:,t], -3),-1)            
    elif op[1] == 'copy':
        raise
    else:
        raise

TMatrices_interp = interpolate.interp1d(freqs_orig*interpfreqfactor, TMatrices, axis=0, kind='linear',fill_value="extrapolate")



# In[4]:

om = np.linspace(np.min(freqs_orig), np.max(freqs_orig),100)
TMatrix0ip = np.reshape(TMatrices_interp(om)[:,0], (len(om), 2*nelem*2*nelem))
f, axa = plt.subplots(2, 2, figsize=(15,15))

print(TMatrices.shape)
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
pdf.savefig(f)


# In[ ]:

kdensity = 66
bz_0 = np.array((0,0,0.,))
bz_K1 = np.array((1.,0,0))*4*np.pi/3/hexside/s3
bz_K2 = np.array((1./2.,s3/2,0))*4*np.pi/3/hexside/s3
bz_M = np.array((3./4, s3/4,0))*4*np.pi/3/hexside/s3
k0Mlist = bz_0 + (bz_M-bz_0) * np.linspace(0,1,kdensity/5)[:,nx]
kMK1list = bz_M + (bz_K1-bz_M) * np.linspace(0,1,kdensity)[:,nx]
kK10list = bz_K1 + (bz_0-bz_K1) * np.linspace(0,1,kdensity)[:,nx]
k0K2list = bz_0 + (bz_K2-bz_0) * np.linspace(0,1,kdensity/5)[:,nx]
kK2Mlist = bz_K2 + (bz_M-bz_K2) * np.linspace(0,1,kdensity/5)[:,nx]
B1 = 2* bz_K1 - bz_K2
B2 = 2* bz_K2 - bz_K1
klist = np.concatenate((k0Mlist,kMK1list,kK10list,k0K2list,kK2Mlist), axis=0)
kxmaplist = np.concatenate((np.array([0]),np.cumsum(np.linalg.norm(np.diff(klist, axis=0), axis=-1))))


# In[ ]:

n2id = np.identity(2*nelem)
n2id.shape = (2,nelem,2,nelem)
extlistlist = list()
leftmatrixlistlist = list()
minsvTElistlist=list()
minsvTMlistlist=list()
nan = float('nan')
omegalist = list()
filecount = 0
for trfile in os.scandir(translations_dir):
    filecount += 1
    try:
        npz = np.load(trfile.path, mmap_mode='r')
        k_0 = npz['precalc_params'][()]['k_hexside'] / hexside
        omega = k_0 * c / math.sqrt(epsilon_b)
        if(omega < 2.4e15 or omega > 2.7e15 ):
            continue
    except:
        print ("Unexpected error, trying to continue with another file:", sys.exc_info()[0])
        continue   
    try:
        tdic = qpms.hexlattice_precalc_AB_loadunwrap(trfile.path, return_points=True)
    except:
        print ("Unexpected error, trying to continue with another file:", sys.exc_info()[0])
        continue
    k_0 = tdic['k_hexside'] / hexside
    omega = k_0 * c / math.sqrt(epsilon_b)
    omegalist.append(omega)
    print(filecount, omega/eV*hbar)
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
        unitcell_envelope = np.exp(-np.sum(unitcell_translations**2,axis=-1)/(2*gaussianSigma**2))
        u2d_envelope = np.exp(-np.sum(u2d_translations**2,axis=-1)/(2*gaussianSigma**2))
        d2u_envelope = np.exp(-np.sum(d2u_translations**2,axis=-1)/(2*gaussianSigma**2))
    
    
    TMatrices_om = TMatrices_interp(omega)

    minsvTElist = np.full((klist.shape[0]),np.nan)
    minsvTMlist = np.full((klist.shape[0]),np.nan)
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
        leftmatrix[0,0,:,1,1,:] = np.tensordot(b_d2u, phases_d2u,axes=(0,-1)) #d2u,M2E
        leftmatrix[0,1,:,1,0,:] = leftmatrix[0,0,:,1,1,:] #d2u, E2M
        leftmatrix[1,0,:,0,1,:] = np.tensordot(b_u2d, phases_u2d,axes=(0,-1)) #u2d,M2E
        leftmatrix[1,1,:,0,0,:] = leftmatrix[1,0,:,0,1,:] #u2d, E2M
        #leftmatrix is now the translation matrix T
        for j in range(2):
            leftmatrix[j] = -np.tensordot(TMatrices_om[j], leftmatrix[j], axes=([-2,-1],[0,1]))
            # at this point, jth row of leftmatrix is that of -MT
            leftmatrix[j,:,:,j,:,:] += n2id
        #now we are done, 1-M

        leftmatrixlist[ki] = leftmatrix


    nnlist = np.logical_not(isNaNlist)
    leftmatrixlist_s = np.reshape(leftmatrixlist,(klist.shape[0], 2*2*nelem,2*2*nelem))[nnlist]
    leftmatrixlist_TE = leftmatrixlist_s[np.ix_(np.arange(leftmatrixlist_s.shape[0]),TEč,TEč)]
    leftmatrixlist_TM = leftmatrixlist_s[np.ix_(np.arange(leftmatrixlist_s.shape[0]),TMč,TMč)]
    minsvTElist[nnlist] = np.amin(np.linalg.svd(leftmatrixlist_TE, compute_uv=False), axis=-1)
    minsvTMlist[nnlist] = np.amin(np.linalg.svd(leftmatrixlist_TM, compute_uv=False), axis=-1)
    minsvTMlistlist.append(minsvTMlist)
    minsvTElistlist.append(minsvTElist) 

    
        


# In[ ]:

minsvTElistarr = np.array(minsvTElistlist)
minsvTMlistarr = np.array(minsvTMlistlist)
omegalist = np.array(omegalist)
omlist = np.broadcast_to(omegalist[:,nx], minsvTElistarr.shape)
kxmlarr = np.broadcast_to(kxmaplist[nx,:], minsvTElistarr.shape)
klist = np.concatenate((k0Mlist,kMK1list,kK10list,k0K2list,kK2Mlist), axis=0)


# In[ ]:

import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as patches
f, ax = plt.subplots(1, figsize=(20,15))
sc = ax.scatter(kxmlarr, omlist/eV*hbar, c = np.sqrt(minsvTMlistarr), s =40, lw=0)
ax.plot(kxmaplist, np.linalg.norm(klist,axis=-1)*cdn/eV*hbar, '-',
       kxmaplist, np.linalg.norm(klist+B1, axis=-1)*cdn/eV*hbar, '-',
       kxmaplist, np.linalg.norm(klist+B2, axis=-1)*cdn/eV*hbar, '-',
       kxmaplist, np.linalg.norm(klist-B2, axis=-1)*cdn/eV*hbar, '-',
       kxmaplist, np.linalg.norm(klist-B1, axis=-1)*cdn/eV*hbar, '-',
       kxmaplist, np.linalg.norm(klist+B2-B1, axis=-1)*cdn/eV*hbar, '-',
        kxmaplist, np.linalg.norm(klist-B2+B1, axis=-1)*cdn/eV*hbar, '-',
        kxmaplist, np.linalg.norm(klist-B2-B1, axis=-1)*cdn/eV*hbar, '-',
        kxmaplist, np.linalg.norm(klist+B2+B1, axis=-1)*cdn/eV*hbar, '-',
        kxmaplist, np.linalg.norm(klist-2*B1, axis=-1)*cdn/eV*hbar, '-',
        kxmaplist, np.linalg.norm(klist-2*B2, axis=-1)*cdn/eV*hbar, '-',
        kxmaplist, np.linalg.norm(klist-2*B2-B1, axis=-1)*cdn/eV*hbar, '-',
        kxmaplist, np.linalg.norm(klist-2*B1-B2, axis=-1)*cdn/eV*hbar, '-',
        kxmaplist, np.linalg.norm(klist-2*B1-2*B2, axis=-1)*cdn/eV*hbar, '-',
#        kxmaplist, np.linalg.norm(klist+2*B2-B1, axis=-1)*cdn, '-',
#        kxmaplist, np.linalg.norm(klist+2*B1-B2, axis=-1)*cdn, '-',
       )
ax.set_xlim([np.min(kxmlarr),np.max(kxmlarr)])
#ax.set_ylim([2.15,2.30])
ax.set_ylim([np.min(omlist/eV*hbar),np.max(omlist/eV*hbar)])
ax.set_xticks([0, kxmaplist[len(k0Mlist)-1], kxmaplist[len(k0Mlist)+len(kMK1list)-1], kxmaplist[len(k0Mlist)+len(kMK1list)+len(kK10list)-1], kxmaplist[len(k0Mlist)+len(kMK1list)+len(kK10list)+len(k0K2list)-1], kxmaplist[len(k0Mlist)+len(kMK1list)+len(kK10list)+len(k0K2list)+len(kK2Mlist)-1]])
ax.set_xticklabels(['Γ', 'M', 'K', 'Γ', 'K\'','M'])
f.colorbar(sc)

pdf.savefig(f)


# In[ ]:

import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as patches
f, ax = plt.subplots(1, figsize=(20,15))
sc = ax.scatter(kxmlarr, omlist/eV*hbar, c = np.sqrt(minsvTElistarr), s =40, lw=0)
ax.plot(kxmaplist, np.linalg.norm(klist,axis=-1)*cdn/eV*hbar, '-',
       kxmaplist, np.linalg.norm(klist+B1, axis=-1)*cdn/eV*hbar, '-',
       kxmaplist, np.linalg.norm(klist+B2, axis=-1)*cdn/eV*hbar, '-',
       kxmaplist, np.linalg.norm(klist-B2, axis=-1)*cdn/eV*hbar, '-',
       kxmaplist, np.linalg.norm(klist-B1, axis=-1)*cdn/eV*hbar, '-',
       kxmaplist, np.linalg.norm(klist+B2-B1, axis=-1)*cdn/eV*hbar, '-',
        kxmaplist, np.linalg.norm(klist-B2+B1, axis=-1)*cdn/eV*hbar, '-',
        kxmaplist, np.linalg.norm(klist-B2-B1, axis=-1)*cdn/eV*hbar, '-',
        kxmaplist, np.linalg.norm(klist+B2+B1, axis=-1)*cdn/eV*hbar, '-',
        kxmaplist, np.linalg.norm(klist-2*B1, axis=-1)*cdn/eV*hbar, '-',
        kxmaplist, np.linalg.norm(klist-2*B2, axis=-1)*cdn/eV*hbar, '-',
        kxmaplist, np.linalg.norm(klist-2*B2-B1, axis=-1)*cdn/eV*hbar, '-',
        kxmaplist, np.linalg.norm(klist-2*B1-B2, axis=-1)*cdn/eV*hbar, '-',
        kxmaplist, np.linalg.norm(klist-2*B1-2*B2, axis=-1)*cdn/eV*hbar, '-',
#        kxmaplist, np.linalg.norm(klist+2*B2-B1, axis=-1)*cdn, '-',
#        kxmaplist, np.linalg.norm(klist+2*B1-B2, axis=-1)*cdn, '-',
       )
ax.set_xlim([np.min(kxmlarr),np.max(kxmlarr)])
#ax.set_ylim([2.15,2.30])
ax.set_ylim([np.min(omlist/eV*hbar),np.max(omlist/eV*hbar)])
ax.set_xticks([0, kxmaplist[len(k0Mlist)-1], kxmaplist[len(k0Mlist)+len(kMK1list)-1], kxmaplist[len(k0Mlist)+len(kMK1list)+len(kK10list)-1], kxmaplist[len(k0Mlist)+len(kMK1list)+len(kK10list)+len(k0K2list)-1], kxmaplist[len(k0Mlist)+len(kMK1list)+len(kK10list)+len(k0K2list)+len(kK2Mlist)-1]])
ax.set_xticklabels(['Γ', 'M', 'K', 'Γ', 'K\'','M'])
f.colorbar(sc)

pdf.savefig(f)


# In[ ]:

pdf.close()


# In[ ]:

unitcell_translations


# In[ ]:



