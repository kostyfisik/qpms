import math
import numpy as np
cimport numpy as np
nx = None

cdef double _s3 = math.sqrt(3)

from scipy.constants import c
from .timetrack import _time_b, _time_e
from .tmatrices import symz_indexarrays
from .hexpoints import hexlattice_get_AB

cpdef hexlattice_zsym_getSVD(int lMax, TMatrices_om, double epsilon_b, double hexside, size_t maxlayer, double omega, klist, gaussianSigma=False, int onlyNmin=0, verbose=False):
    cdef np.ndarray[np.npy_double, ndim = 2] klist_c = klist
    btime = _time_b(verbose)
    cdef size_t nelem = lMax * (lMax + 2)
    _n2id = np.identity(2*nelem)
    _n2id.shape = (2,nelem,2,nelem)
    cdef np.ndarray[np.npy_double, ndim = 4] n2id = _n2id
    cdef double nan = float('nan')
    k_0 = omega * math.sqrt(epsilon_b) / c
    tdic = hexlattice_get_AB(lMax,k_0*hexside,maxlayer)
    cdef np.ndarray[np.npy_cdouble, ndim = 3] a_self = tdic['a_self'][:,:nelem,:nelem]
    cdef np.ndarray[np.npy_cdouble, ndim = 3] b_self = tdic['b_self'][:,:nelem,:nelem]
    cdef np.ndarray[np.npy_cdouble, ndim = 3] a_u2d = tdic['a_u2d'][:,:nelem,:nelem]
    cdef np.ndarray[np.npy_cdouble, ndim = 3] b_u2d = tdic['b_u2d'][:,:nelem,:nelem]
    cdef np.ndarray[np.npy_cdouble, ndim = 3] a_d2u = tdic['a_d2u'][:,:nelem,:nelem]
    cdef np.ndarray[np.npy_cdouble, ndim = 3] b_d2u = tdic['b_d2u'][:,:nelem,:nelem]
    cdef np.ndarray[np.npy_double, ndim = 2] unitcell_translations = tdic['self_tr']*hexside*_s3
    cdef np.ndarray[np.npy_double, ndim = 2] u2d_translations = tdic['u2d_tr']*hexside*_s3
    cdef np.ndarray[np.npy_double, ndim = 2] d2u_translations = tdic['d2u_tr']*hexside*_s3

    cdef np.ndarray[np.npy_double, ndim = 1] unitcell_envelope, u2d_envelope, d2u_envelope
    if gaussianSigma:
        sbtime = _time_b(verbose, step='Calculating gaussian envelope')
        unitcell_envelope = np.exp(-np.sum(tdic['self_tr']**2,axis=-1)/(2*gaussianSigma**2))
        u2d_envelope = np.exp(-np.sum(tdic['u2d_tr']**2,axis=-1)/(2*gaussianSigma**2))
        d2u_envelope = np.exp(-np.sum(tdic['d2u_tr']**2,axis=-1)/(2*gaussianSigma**2))
        _time_e(sbtime, verbose, step='Calculating gaussian envelope')

    cdef np.ndarray[np.npy_cdouble, ndim = 3] svUfullTElist, svVfullTElist, svUfullTMlist, svVfullTMlist
    cdef np.ndarray[np.npy_cdouble, ndim = 2] svSfullTElist, svSfullTMlist
    cdef np.ndarray[np.npy_double, ndim = 2] minsvTElist, minsvTMlist
    
    #TMatrices_om = TMatrices_interp(omega)
    if(not onlyNmin):
        svUfullTElist = np.full((klist_c.shape[0], 2*nelem, 2*nelem), np.nan, dtype=complex)
        svVfullTElist = np.full((klist_c.shape[0], 2*nelem, 2*nelem), np.nan, dtype=complex)
        svSfullTElist = np.full((klist_c.shape[0], 2*nelem), np.nan, dtype=complex)
        svUfullTMlist = np.full((klist_c.shape[0], 2*nelem, 2*nelem), np.nan, dtype=complex)
        svVfullTMlist = np.full((klist_c.shape[0], 2*nelem, 2*nelem), np.nan, dtype=complex)
        svSfullTMlist = np.full((klist_c.shape[0], 2*nelem), np.nan, dtype=complex)
    else:
        minsvTElist = np.full((klist_c.shape[0], onlyNmin),np.nan)
        minsvTMlist = np.full((klist_c.shape[0], onlyNmin),np.nan)

    cdef np.ndarray[np.npy_cdouble] leftmatrixlist = np.full((klist_c.shape[0],2,2,nelem,2,2,nelem),np.nan,dtype=complex)
    cdef np.ndarray[np.npy_bool, ndim=1] isNaNlist = np.zeros((klist_c.shape[0]), dtype=bool)

    sbtime = _time_b(verbose, step='Initialization of matrices for SVD for a given list of k\'s')
    # sem nějaká rozumná smyčka

    # declarations for the ki loop:
    cdef size_t ki
    cdef np.ndarray[np.npy_cdouble, ndim = 1] phases_self
    cdef np.ndarray[np.npy_cdouble, ndim = 1] phases_u2d 
    cdef np.ndarray[np.npy_cdouble, ndim = 1] phases_d2u 
    cdef np.ndarray[np.npy_cdouble, ndim=6] leftmatrix 
    cdef np.ndarray[np.npy_double, ndim=1] k
    cdef int j
    
    for ki in range(klist_c.shape[0]):
        k = klist_c[ki]
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
        #       0:[u,E<--u,E  ]
        #       1:[d,M<--d,M  ]
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
    leftmatrixlist_s = np.reshape(leftmatrixlist,(klist_c.shape[0], 2*2*nelem,2*2*nelem))[nnlist]
    TEc, TMc = symz_indexarrays(lMax, 2)
    leftmatrixlist_TE = leftmatrixlist_s[np.ix_(np.arange(leftmatrixlist_s.shape[0]),TEc,TEc)]
    leftmatrixlist_TM = leftmatrixlist_s[np.ix_(np.arange(leftmatrixlist_s.shape[0]),TMc,TMc)]
    _time_e(sbtime, verbose, step='Initializing matrices for SVD for a given list of k\'s')

    sbtime = _time_b(verbose, step='Calculating SVDs for a given list of k\'s.')
    if(not onlyNmin):
        svUfullTElist[nnlist], svSfullTElist[nnlist], svVfullTElist[nnlist] = np.linalg.svd(leftmatrixlist_TE, compute_uv=True)
        svUfullTMlist[nnlist], svSfullTMlist[nnlist], svVfullTMlist[nnlist] = np.linalg.svd(leftmatrixlist_TM, compute_uv=True)
        _time_e(sbtime, verbose, step='Calculating SVDs for a given list of k\'s.')
        return ((svUfullTElist, svSfullTElist, svVfullTElist), (svUfullTMlist, svSfullTMlist, svVfullTMlist))
    else:
        minsvTElist[nnlist] = np.linalg.svd(leftmatrixlist_TE, compute_uv=False)[...,-onlyNmin:]
        minsvTMlist[nnlist] = np.linalg.svd(leftmatrixlist_TM, compute_uv=False)[...,-onlyNmin:]
        _time_e(sbtime, verbose, step='Calculating SVDs for a given list of k\'s.')
        return (minsvTElist, minsvTMlist)
