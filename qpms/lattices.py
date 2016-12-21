'''
Object oriented approach for the classical multiple scattering problem.
'''
import numpy as np
import time
import scipy
import sys
from qpms_c import get_mn_y # TODO be explicit about what is imported
from .qpms_p import cart2sph, nelem2lMax, Ã, B̃ # TODO be explicit about what is imported

def _time_b(active = True, name = None, step = None):
    '''
    Auxiliary function for keeping track of elapsed time.
    Returns current time (to be used by _time_e).
    '''
    now = time.time()
    if active:
        if not name:
            name = sys._getframe(1).f_code.co_name
        if step:
            print('%.4f: %s in function %s started.' % (now, step, name), file = sys.stderr)
        else:
            print('%.4f: Function %s started.' % (now, name), file=sys.stderr)
        sys.stderr.flush()
    return now

def _time_e(start_time, active = True, name = None, step = None):
    now = time.time()
    if active:
        if not name:
            name = sys._getframe(1).f_code.co_name
        if step:
            print('%.4f: %s in function %s finished (elapsed %.2f s).' 
                    % (now, step, name, now - start_time), file = sys.stderr)
        else:
            print('%.4f: Function %s finished (elapsed %.2f s).' 
                    % (now, name, now - start_time), file = sys.stderr)
        sys.stderr.flush()

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
        self.positions = positions
        self.interaction_matrix = None
        self.N = positions.shape[0]
        self.k_0 = k_0
        self.lMax = lMax if lMax else nelem2lMax(TMatrices.shape[-1])
        nelem = lMax * (lMax + 2) #!
        self.nelem = nelem #!
        self.prepared = False
        self.TMatrices = np.broadcast_to(TMatrices, (self.N,2,nelem,2,nelem))

    def prepare(self, keep_interaction_matrix = False, verbose=False):
        btime = _time_b(verbose)
        if not self.prepared:
            if not self.interaction_matrix:
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
        self.prepare(verbose=verbose)
        pq_0 = np.broadcast_to(pq_0, (self.N,2,self.nelem))
        MP_0 = np.empty((N,2,nelem),dtype=np.complex_)
        for j in range(self.N):
            MP_0[j] = np.tensordot(self.TMatrices[j], pq_0[j],axes=[-2,-1],[-2,-1])
        MP_0.shape = (N*2*self.nelem,)
        solvebtime = _time_b(verbose,step='Solving the linear equation')
        ab = scipy.linalg.lu_solve(self.lupiv, MP_0)
        _time_e(solvebtime, verbose, step='Solving the linear equation')
        ab.shape = (N,2,nelem)
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

class Scattering_2D_lattice(Scattering):
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
        self.TMatrices = np.broadcast_to(TMatrices, (self.N,2,nelem,2,nelem))   def __init__(self):
