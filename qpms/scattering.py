'''
Object oriented approach for the classical multiple scattering problem.
'''

__TODO__ = '''
- Implement per-scatterer lMax
  - This means that Scattering.TMatrices either can not be a single array with a fixed
    (N, 2, nelem, 2, nelem) shape but rather list of (2, nelem, 2, nelem) with nelem varying
    per particle or some of its elements have to be unused. Anyways, there has to be some kind of
    list with the lMaxes.
'''

import numpy as np
nx = np.newaxis
import time
import scipy
import sys
import warnings
import math
from qpms_c import get_mn_y, trans_calculator # TODO be explicit about what is imported
from .qpms_p import cart2sph, nelem2lMax # TODO be explicit about what is imported
from .timetrack import _time_b, _time_e

class Scattering(object):
    '''

    This is the most general class for a system of scatterers
    in a non-lossy homogeneous background
    to be solved with the multiple_scattering method. The scatterers,
    as long as they comply with the disjoint circumscribed sphere
    hypothesis, can each have any position in the 3D space and
    any T-matrix.

    Note that this object describes the scattering problem only for
    a single given frequency, as the T-matrices and wavelenght
    otherwise differ and all the computationally demanding
    parts have to be done for each frequency. However,
    the object can be recycled for many incident field shapes
    at the given frequency.

    Attributes should be perhaps later redefined to be read-only
    (or make descriptors for them).

    Args:
        positions: (N,3)-shaped real array
        TMatrices: (N,2,nelem,2,nelem)-shaped array
        k_0 (float): Wave number for the space between scatterers.

    Attributes:
        positions:
        TMatrices:
        k_0 (float): Wave number for the space between scatterers.
        lMax (int): Absolute maximum l for all scatterers. Depending on implementation,
            lMax can be smaller for some individual scatterers in certain subclasses.
            FIXME: here it is still implemented as constant lMax for all sites, see #!
        prepared (bool): Keeps information whether the interaction matrix has
            already been built and factorized.
        

    '''

    def __init__(self, positions, TMatrices, k_0, lMax = None, verbose=False, J_scat=3):
        self.J_scat = J_scat
        self.positions = positions.reshape((-1, positions.shape[-1]))
        self.interaction_matrix = None
        self.N = self.positions.shape[0]
        self.k_0 = k_0
        self.lMax = lMax if lMax else nelem2lMax(TMatrices.shape[-1])
        self.tc = trans_calculator(self.lMax)
        nelem = self.lMax * (self.lMax + 2) #!
        self.nelem = nelem #!
        self.prepared = False
        self.TMatrices = np.broadcast_to(TMatrices, (self.N,2,nelem,2,nelem))
        if np.isnan(np.min(TMatrices)):
            warnings.warn("Some TMatrices contain NaNs. Expect invalid results")
        if np.isnan(np.min(positions)):
            warnings.warn("positions contain NaNs. Expect invalid results")
        if math.isnan(k_0):
            warnings.warn("k_0 is NaN. Expect invalid results")



    def prepare(self, keep_interaction_matrix = False, verbose=False):
        btime = _time_b(verbose)
        if not self.prepared:
            if self.interaction_matrix is None:
                self.build_interaction_matrix(verbose=verbose)
            self.lupiv = scipy.linalg.lu_factor(self.interaction_matrix,overwrite_a = not keep_interaction_matrix)
            if not keep_interaction_matrix:
                self.interaction_matrix = None
            self.prepared = True
        _time_e(btime, verbose)

    def build_interaction_matrix(self,verbose = False):
        btime = _time_b(verbose)
        N = self.N
        my, ny = get_mn_y(self.lMax)
        nelem = len(my)
        leftmatrix = np.zeros((N,2,nelem,N,2,nelem), dtype=complex)
        sbtime = _time_b(verbose, step = 'Calculating interparticle translation coefficients')
        """
        for i in range(N):
            for j in range(N):
                for yi in range(nelem):
                    for yj in range(nelem):
                        if(i != j):
                            d_i2j = cart2sph(self.positions[j]-self.positions[i])
                            a = Ã(my[yj],ny[yj],my[yi],ny[yi],kdlj=d_i2j[0]*self.k_0,θlj=d_i2j[1],φlj=d_i2j[2],r_ge_d=False,J=self.J_scat)
                            b = B̃(my[yj],ny[yj],my[yi],ny[yi],kdlj=d_i2j[0]*self.k_0,θlj=d_i2j[1],φlj=d_i2j[2],r_ge_d=False,J=self.J_scat)
                            leftmatrix[j,0,yj,i,0,yi] = a
                            leftmatrix[j,1,yj,i,1,yi] = a
                            leftmatrix[j,0,yj,i,1,yi] = b
                            leftmatrix[j,1,yj,i,0,yi] = b
        """
        kdji = cart2sph(self.positions[:,nx,:] - self.positions[nx,:,:])
        kdji[:,:,0] *= self.k_0
        # get_AB array structure: [j,yj,i,yi]
        a, b = self.tc.get_AB(my[nx,:,nx,nx],ny[nx,:,nx,nx],my[nx,nx,nx,:],ny[nx,nx,nx,:],
                (kdji[:,:,0])[:,nx,:,nx], (kdji[:,:,1])[:,nx,:,nx], (kdji[:,:,2])[:,nx,:,nx],
                False,self.J_scat)
        mask = np.broadcast_to(np.eye(N,dtype=bool)[:,nx,:,nx],(N,nelem,N,nelem))
        a[mask] = 0 # no self-translations
        b[mask] = 0 
        leftmatrix[:,0,:,:,0,:] = a
        leftmatrix[:,1,:,:,1,:] = a
        leftmatrix[:,0,:,:,1,:] = b
        leftmatrix[:,1,:,:,0,:] = b
        _time_e(sbtime, verbose, step = 'Calculating interparticle translation coefficients')
        # at this point, leftmatrix is the translation matrix
        n2id = np.identity(2*nelem)
        n2id.shape = (2,nelem,2,nelem)
        for j in range(N):
            leftmatrix[j] = - np.tensordot(self.TMatrices[j],leftmatrix[j],axes=([-2,-1],[0,1]))
            # at this point, jth row of leftmatrix is that of -MT
            leftmatrix[j,:,:,j,:,:] += n2id
            # now we are done, 1-MT
        leftmatrix.shape=(N*2*nelem,N*2*nelem)
        self.interaction_matrix = leftmatrix
        _time_e(btime, verbose)

    def scatter(self, pq_0, verbose = False):
        '''
        pq_0 is (N, nelem, 2)-shaped array 
        '''
        btime = _time_b(verbose)
        if math.isnan(np.min(pq_0)):
            warnings.warn("The incident field expansion contains NaNs. Expect invalid results.")
        self.prepare(verbose=verbose)
        pq_0 = np.broadcast_to(pq_0, (self.N,2,self.nelem))
        MP_0 = np.empty((self.N,2,self.nelem),dtype=np.complex_)
        for j in range(self.N):
            MP_0[j] = np.tensordot(self.TMatrices[j], pq_0[j],axes=([-2,-1],[-2,-1]))
        MP_0.shape = (self.N*2*self.nelem,)
        solvebtime = _time_b(verbose,step='Solving the linear equation')
        ab = scipy.linalg.lu_solve(self.lupiv, MP_0)
        if math.isnan(np.min(ab)):
            warnings.warn("Got NaN in the scattering result. Damn.")
            raise
        _time_e(solvebtime, verbose, step='Solving the linear equation')
        ab.shape = (self.N,2,self.nelem)
        _time_e(btime, verbose)
        return ab
        
    def scatter_constmultipole(self, pq_0_c, verbose = False):
        btime = _time_b(verbose)
        N = self.N
        self.prepare(verbose=verbose)
        nelem = self.nelem
        if(pq_0_c ==1):
            pq_0_c = np.full((2,nelem),1)
        ab = np.empty((2,nelem,N*2*nelem), dtype=complex)
        for N_or_M in range(2):
            for yy in range(nelem):
                pq_0 = np.zeros((2,nelem),dtype=np.complex_)
                pq_0[N_or_M,yy] = pq_0_c[N_or_M,yy]
                pq_0 = np.broadcast_to(pq_0, (N,2,nelem))
                MP_0 = np.empty((N,2,nelem),dtype=np.complex_)
                for j in range(N):
                    MP_0[j] = np.tensordot(self.TMatrices[j], pq_0[j],axes=([-2,-1],[-2,-1]))
                MP_0.shape = (N*2*nelem,)
                ab[N_or_M,yy] = scipy.linalg.lu_solve(self.lupiv,MP_0)
        ab.shape = (2,nelem,N,2,nelem)
        _time_e(btime, verbose)
        return ab

class LatticeScattering(Scattering):
    def __init__(self, lattice_spec, k_0, zSym = False):


"""
class Scattering_2D_lattice_rectcells(Scattering):
    def __init__(self, rectcell_dims, rectcell_elem_positions, cellspec, k_0, rectcell_TMatrices = None, TMatrices = None, lMax = None, verbose=False, J_scat=3):
        '''
        cellspec: dvojice ve tvaru (seznam_zaplněnosti, seznam_pozic)
        '''
        if (rectcell_TMatrices is None) == (TMatrices is None):
            raise ValueError('Either rectcell_TMatrices or TMatrices has to be given')
        ###self.positions = ZDE JSEM SKONČIL
        self.J_scat = J_scat
        self.positions = positions
        self.interaction_matrix = None
        self.N = positions.shape[0]
        self.k_0 = k_0
        self.lMax = lMax if lMax else nelem2lMax(TMatrices.shape[-1])
        nelem = lMax * (lMax + 2) #!
        self.nelem = nelem #!
        self.prepared = False
        self.TMatrices = np.broadcast_to(TMatrices, (self.N,2,nelem,2,nelem))   
"""

class Scattering_2D_zsym(Scattering):
    def __init__(self, positions, TMatrices, k_0, lMax = None, verbose=False, J_scat=3):
        Scattering.__init__(self, positions, TMatrices, k_0, lMax, verbose, J_scat)
        #TODO some checks on TMatrices symmetry
        self.TE_yz = np.arange(self.nelem)
        self.TM_yz = self.TE_yz
        self.my, self.ny = get_mn_y(self.lMax)
        self.TE_NMz = (self.my + self.ny) % 2
        self.TM_NMz = 1 - self.TE_NMz
        self.tc = trans_calculator(self.lMax)
        # TODO možnost zadávat T-matice rovnou ve zhuštěné podobě
        TMatrices_TE = TMatrices[...,self.TE_NMz[:,nx],self.TE_yz[:,nx],self.TE_NMz[nx,:],self.TE_yz[nx,:]]
        TMatrices_TM = TMatrices[...,self.TM_NMz[:,nx],self.TM_yz[:,nx],self.TM_NMz[nx,:],self.TM_yz[nx,:]]
        self.TMatrices_TE = np.broadcast_to(TMatrices_TE, (self.N, self.nelem, self.nelem))
        self.TMatrices_TM = np.broadcast_to(TMatrices_TM, (self.N, self.nelem, self.nelem))
        self.prepared_TE = False
        self.prepared_TM = False
        self.interaction_matrix_TE = None
        self.interaction_matrix_TM= None

    def prepare_partial(self, TE_or_TM, keep_interaction_matrix = False, verbose=False):
        '''
        TE is 0, TM is 1.
        '''
        btime = _time_b(verbose)
        if (TE_or_TM == 0): #TE
            if not self.prepared_TE:
                if self.interaction_matrix_TE is None:
                    self.build_interaction_matrix(0, verbose)
                sbtime = _time_b(verbose, step = 'Calculating LU decomposition of the interaction matrix, TE part')
                self.lupiv_TE = scipy.linalg.lu_factor(self.interaction_matrix_TE, overwrite_a = not keep_interaction_matrix)
                _time_e(sbtime, verbose, step = 'Calculating LU decomposition of the interaction matrix, TE part')
                if(np.isnan(np.min(self.lupiv_TE[0])) or np.isnan(np.min(self.lupiv_TE[1]))):
                    warnings.warn("LU factorisation of interaction matrix contains NaNs. Expect invalid results.")
                self.prepared_TE = True
        if (TE_or_TM == 1): #TM
            if not self.prepared_TM:
                if self.interaction_matrix_TM is None:
                    self.build_interaction_matrix(1, verbose)
                sbtime = _time_b(verbose, step = 'Calculating LU decomposition of the interaction matrix, TM part')
                self.lupiv_TM = scipy.linalg.lu_factor(self.interaction_matrix_TM, overwrite_a = not keep_interaction_matrix)
                _time_e(sbtime, verbose, step = 'Calculating LU decomposition of the interaction matrix, TM part')
                if(np.isnan(np.min(self.lupiv_TM[0])) or np.isnan(np.min(self.lupiv_TM[1]))):
                    warnings.warn("LU factorisation of interaction matrix contains NaNs. Expect invalid results.")
                self.prepared_TM = True
        _time_e(btime, verbose)

    def prepare(self, keep_interaction_matrix = False, verbose=False):
        btime = _time_b(verbose)
        if not self.prepared:
            self.prepare_partial(0, keep_interaction_matrix, verbose)
            self.prepare_partial(1, keep_interaction_matrix, verbose)
            self.prepared = True
        _time_e(btime, verbose)

    def build_interaction_matrix(self,TE_or_TM = None, verbose = False):
        #None means both
        btime = _time_b(verbose)
        N = self.N
        my, ny = get_mn_y(self.lMax)
        nelem = len(my)
        idm = np.identity(nelem)
        if (TE_or_TM == 0):
            EoMl = (0,)
        elif (TE_or_TM == 1):
            EoMl = (1,)
        elif (TE_or_TM is None):
            EoMl = (0,1)
        sbtime = _time_b(verbose, step = 'Calculating interparticle translation coefficients')
        kdji = cart2sph(self.positions[:,nx,:] - self.positions[nx,:,:], allow2d=True)
        kdji[:,:,0] *= self.k_0
        # get_AB array structure: [j,yj,i,yi]
        # FIXME I could save some memory by calculating only half of these coefficients
        a, b = self.tc.get_AB(my[nx,:,nx,nx],ny[nx,:,nx,nx],my[nx,nx,nx,:],ny[nx,nx,nx,:],
                (kdji[:,:,0])[:,nx,:,nx], (kdji[:,:,1])[:,nx,:,nx], (kdji[:,:,2])[:,nx,:,nx],
                False,self.J_scat)
        mask = np.broadcast_to(np.eye(N,dtype=bool)[:,nx,:,nx],(N,nelem,N,nelem))
        a[mask] = 0 # no self-translations
        b[mask] = 0
        if np.isnan(np.min(a)) or np.isnan(np.min(b)):
            warnings.warn("Some of the translation coefficients is a NaN. Expect invalid results.")
        _time_e(sbtime, verbose, step = 'Calculating interparticle translation coefficients')
        for EoM in EoMl:
            leftmatrix = np.zeros((N,nelem,N,nelem), dtype=complex)
            y = np.arange(nelem)
            yi = y[nx,nx,nx,:]
            yj = y[nx,:,nx,nx]
            mask = np.broadcast_to((((yi - yj) % 2) == 0),(N,nelem,N,nelem))
            leftmatrix[mask] = a[mask]
            mask = np.broadcast_to((((yi - yj) % 2) != 0),(N,nelem,N,nelem))
            leftmatrix[mask] = b[mask]
            """ # we use to calculate the AB coefficients here
            for i in range(N):
                for j in range(i):
                    for yi in range(nelem):
                        for yj in range(nelem):
                            d_i2j = cart2sph(self.positions[j]-self.positions[i])
                            if ((yi - yj) % 2) == 0:
                                tr = Ã(my[yj],ny[yj],my[yi],ny[yi],kdlj=d_i2j[0]*self.k_0,θlj=d_i2j[1],φlj=d_i2j[2],r_ge_d=False,J=self.J_scat)
                            else:
                                tr = B̃(my[yj],ny[yj],my[yi],ny[yi],kdlj=d_i2j[0]*self.k_0,θlj=d_i2j[1],φlj=d_i2j[2],r_ge_d=False,J=self.J_scat)
                            leftmatrix[j,yj,i,yi] = tr
                            leftmatrix[i,yi,j,yj] = tr if (0 == (my[yj]+my[yi]) % 2) else -tr
            _time_e(sbtime, verbose, step = 'Calculating interparticle translation coefficients, T%s part' % ('M' if EoM else 'E'))
            """
            for j in range(N):
                leftmatrix[j] = - np.tensordot(self.TMatrices_TM[j] if EoM else self.TMatrices_TE[j],leftmatrix[j],
                        axes = ([-1],[0]))
                leftmatrix[j,:,j,:] += idm
            leftmatrix.shape = (self.N*self.nelem, self.N*self.nelem)
            if np.isnan(np.min(leftmatrix)):
                warnings.warn("Interaction matrix contains some NaNs. Expect invalid results.")
            if EoM == 0:
                self.interaction_matrix_TE = leftmatrix
            if EoM == 1:
                self.interaction_matrix_TM = leftmatrix
        a = None
        b = None
        _time_e(btime, verbose)

    def scatter_partial(self, TE_or_TM, pq_0, verbose = False):
        '''
        pq_0 is (N, nelem)-shaped array
        '''
        btime = _time_b(verbose)
        self.prepare_partial(TE_or_TM, verbose = verbose)
        pq_0 = np.broadcast_to(pq_0, (self.N, self.nelem))
        MP_0 = np.empty((self.N,self.nelem),dtype=np.complex_)
        for j in range(self.N):
            if TE_or_TM: #TM
                MP_0[j] = np.tensordot(self.TMatrices_TM[j], pq_0[j], axes=([-1],[-1]))
            else: #TE
                MP_0[j] = np.tensordot(self.TMatrices_TE[j], pq_0[j], axes=([-1],[-1]))
        MP_0.shape = (self.N*self.nelem,)
        solvebtime = _time_b(verbose,step='Solving the linear equation')
        ab = scipy.linalg.lu_solve(self.lupiv_TM if TE_or_TM else self.lupiv_TE, MP_0)
        _time_e(solvebtime, verbose, step='Solving the linear equation')
        ab.shape = (self.N, self.nelem)
        _time_e(btime,verbose)
        return ab

    def scatter(self, pq_0, verbose = False):
        '''
        FI7ME
        pq_0 is (N, nelem, 2)-shaped array 
        '''
        btime = _time_b(verbose)
        raise Exception('Not yet implemented')
        self.prepare(verbose=verbose)
        pq_0 = np.broadcast_to(pq_0, (self.N,2,self.nelem))
        MP_0 = np.empty((N,2,nelem),dtype=np.complex_)
        for j in range(self.N):
            MP_0[j] = np.tensordot(self.TMatrices[j], pq_0[j],axes=([-2,-1],[-2,-1]))
        MP_0.shape = (N*2*self.nelem,)
        solvebtime = _time_b(verbose,step='Solving the linear equation')
        ab = scipy.linalg.lu_solve(self.lupiv, MP_0)
        _time_e(solvebtime, verbose, step='Solving the linear equation')
        ab.shape = (N,2,nelem)
        _time_e(btime, verbose)
        return ab

    def forget_matrices(self):
        '''
        Free interaction matrices and set the respective flags 
        (useful when memory is a bottleneck).
        '''
        self.interaction_matrix_TE = None
        self.interaction_matrix_TM = None
        self.lupiv_TE = None
        self.lupiv_TM = None
        self.prepared_TE = False
        self.prepared_TM = False
        self.prepared = False


