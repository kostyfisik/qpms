"""@package qpms_c
            self.s.norm = <qpms_normalisation_t>(QPMS_NORMALISATION_NORM_POWER | QPMS_NORMALISATION_CSPHASE)
Cythonized parts of QPMS; mostly wrappers over the C data structures
to make them available in Python.
"""

# Cythonized parts of QPMS here
# -----------------------------

import numpy as np
import cmath
from qpms_cdefs cimport *
from cyquaternions cimport *
#from cyquaternions import *
from cybspec cimport *
#from cybspec import *
from cycommon import *
from cycommon cimport make_c_string
cimport cython
import enum
import warnings
import os
from libc.stdlib cimport malloc, free, calloc, abort

cdef class TMatrixInterpolator:
    '''
    Wrapper over the qpms_tmatrix_interpolator_t structure.
    '''
    #cdef readonly np.ndarray m # Numpy array holding the matrix data
    cdef readonly BaseSpec spec # Here we hold the base spec for the correct reference counting; TODO check if it gets copied
    cdef qpms_tmatrix_t *tmatrices_array
    cdef cdouble *tmdata
    cdef double *freqs
    cdef double *freqs_su
    cdef size_t nfreqs
    cdef qpms_tmatrix_interpolator_t *interp

    def __cinit__(self, filename, BaseSpec bspec,  *args, **kwargs):
        '''Creates a T-matrix interpolator object from a scuff-tmatrix output'''
        global qpms_load_scuff_tmatrix_crash_on_failure
        qpms_load_scuff_tmatrix_crash_on_failure = False
        self.spec = bspec
        cdef char * cpath = make_c_string(filename)
        retval = qpms_load_scuff_tmatrix(cpath, self.spec.rawpointer(),
                &(self.nfreqs), &(self.freqs), &(self.freqs_su),
                &(self.tmatrices_array), &(self.tmdata))
        if (retval != QPMS_SUCCESS):
            raise IOError("Could not read T-matrix from %s: %s" % (filename, os.strerror(retval)))
        if 'symmetrise' in kwargs:
            sym = kwargs['symmetrise']
            if isinstance(sym, FinitePointGroup):
                if QPMS_SUCCESS != qpms_symmetrise_tmdata_finite_group(
                        self.tmdata, self.nfreqs, self.spec.rawpointer(),
                        (<FinitePointGroup?>sym).rawpointer()):
                    raise Exception("This should not happen.")
                atol = kwargs['atol'] if 'atol' in kwargs else 1e-16
                qpms_czero_roundoff_clean(self.tmdata, self.nfreqs * len(bspec)**2, atol)
            else:
                warnings.warn('symmetrise argument type not supported; ignoring.')
        self.interp = qpms_tmatrix_interpolator_create(self.nfreqs,
                self.freqs, self.tmatrices_array, gsl_interp_cspline)
        if not self.interp: raise Exception("Unexpected NULL at interpolator creation.")
    def __call__(self, double freq):
        '''Returns a TMatrix instance, corresponding to a given frequency.'''
        if freq < self.freqs[0] or freq > self.freqs[self.nfreqs-1]:# FIXME here I assume that the input is already sorted
            raise ValueError("input frequency %g is outside the interpolator domain (%g, %g)"
                    % (freq, self.freqs[0], self.freqs[self.nfreqs-1]))
        # This is a bit stupid, I should rethink the CTMatrix constuctors
        cdef qpms_tmatrix_t *t = qpms_tmatrix_interpolator_eval(self.interp, freq)
        cdef CTMatrix res = CTMatrix(self.spec, <cdouble[:len(self.spec),:len(self.spec)]>(t[0].m))
        qpms_tmatrix_free(t)
        return res
    def __dealloc__(self):
        qpms_tmatrix_interpolator_free(self.interp)
        free(self.tmatrices_array)
        free(self.tmdata)
        free(self.freqs_su)
        free(self.freqs)
    property freq_interval:
        def __get__(self):
            return [self.freqs[0], self.freqs[self.nfreqs-1]]

cdef class CTMatrix: # N.B. there is another type called TMatrix in tmatrices.py!
    '''
    Wrapper over the C qpms_tmatrix_t stucture. 
    '''
    cdef readonly np.ndarray m # Numpy array holding the matrix data
    cdef readonly BaseSpec spec # Here we hold the base spec for the correct reference counting; TODO check if it gets copied
    cdef qpms_tmatrix_t t

    def __cinit__(CTMatrix self, BaseSpec spec, matrix):
        self.spec = spec
        self.t.spec = self.spec.rawpointer();
        if (matrix is None) or not np.any(matrix):
            self.m = np.zeros((len(spec),len(spec)), dtype=complex, order='C')
        else:
            # The following will raise an exception if shape is wrong
            self.m = np.array(matrix, dtype=complex, copy=True, order='C').reshape((len(spec), len(spec)))
        #self.m.setflags(write=False) # checkme
        cdef cdouble[:,:] m_memview = self.m
        self.t.m = &(m_memview[0,0])
        self.t.owns_m = False # Memory in self.t.m is "owned" by self.m, not by self.t...

    cdef qpms_tmatrix_t *rawpointer(CTMatrix self):
        '''Pointer to the qpms_tmatrix_t structure.
        Don't forget to reference the BaseSpec object itself when storing the pointer anywhere!!!
        '''
        return &(self.t)
    property rawpointer:
        def __get__(self):
            return <uintptr_t> &(self.t)

    # Transparent access to the T-matrix elements.
    def __getitem__(self, key):
        return self.m[key]
    def __setitem__(self, key, value):
        self.m[key] = value

    def as_ndarray(CTMatrix self):
        ''' Returns a copy of the T-matrix as a numpy array.'''
        # Maybe not totally needed after all, as np.array(T[...]) should be equivalent and not longer
        return np.array(self.m, copy=True)

    def spherical_fill(CTMatrix self, double radius, cdouble k_int,
            cdouble k_ext, cdouble mu_int = 1, cdouble mu_ext = 1):
        ''' Replaces the contents of the T-matrix with those of a spherical particle.'''
        qpms_tmatrix_spherical_fill(&self.t, radius, k_int, k_ext, mu_int, mu_ext)

    def spherical_perm_fill(CTMatrix self, double radius, double freq, cdouble epsilon_int,
            cdouble epsilon_ext):
        '''Replaces the contenst of the T-matrix with those of a spherical particle.'''
        qpms_tmatrix_spherical_mu0_fill(&self.t, radius, freq, epsilon_int, epsilon_ext)
        
    @staticmethod
    def spherical(BaseSpec spec, double radius, cdouble k_int, cdouble k_ext, 
            cdouble mu_int = 1, cdouble mu_ext = 1):
        ''' Creates a T-matrix of a spherical nanoparticle. '''
        tm = CTMatrix(spec, 0)
        tm.spherical_fill(radius, k_int, k_ext, mu_int, mu_ext)
        return tm
    
    @staticmethod
    def spherical_perm(BaseSpec spec, double radius, double freq, cdouble epsilon_int, cdouble epsilon_ext):
        '''Creates a T-matrix of a spherical nanoparticle.'''
        tm = CTMatrix(spec, 0)
        tm.spherical_perm_fill(radius, freq, epsilon_int, epsilon_ext)
        return tm

cdef class PointGroup:
    cdef readonly qpms_pointgroup_t G

    def __init__(self, cls, qpms_gmi_t n = 0, IRot3 orientation = IRot3()):
        cls = PointGroupClass(cls)
        self.G.c = cls
        if n <= 0 and qpms_pg_is_finite_axial(cls):
            raise ValueError("For finite axial groups, n argument must be positive")
        self.G.n = n
        self.G.orientation = orientation.qd

    def __len__(self):
        return qpms_pg_order(self.G.c, self.G.n);

    def __le__(PointGroup self, PointGroup other):
        if qpms_pg_is_subgroup(self.G, other.G):
            return True
        else:
            return False
    def __ge__(PointGroup self, PointGroup other):
        if qpms_pg_is_subgroup(other.G, self.G):
            return True
        else:
            return False
    def __lt__(PointGroup self, PointGroup other):
        return qpms_pg_is_subgroup(self.G, other.G) and not qpms_pg_is_subgroup(other.G, self.G)
    def __eq__(PointGroup self, PointGroup other):
        return qpms_pg_is_subgroup(self.G, other.G) and qpms_pg_is_subgroup(other.G, self.G)
    def __gt__(PointGroup self, PointGroup other):
        return not qpms_pg_is_subgroup(self.G, other.G) and qpms_pg_is_subgroup(other.G, self.G)

    def elems(self):
        els = list()
        cdef qpms_irot3_t *arr
        arr = qpms_pg_elems(NULL, self.G)
        cdef IRot3 q
        for i in range(len(self)):
            q = IRot3()
            q.cset(arr[i])
            els.append(q)
        free(arr)
        return els

cdef class FinitePointGroup:
    '''
    Wrapper over the qpms_finite_group_t structure.

    TODO more functionality to make it better usable in Python
    (group element class at least)
    '''
    cdef readonly bint owns_data
    cdef qpms_finite_group_t *G

    def __cinit__(self, info):
        '''Constructs a FinitePointGroup from PointGroupInfo'''
        # TODO maybe I might use a try..finally statement to avoid leaks
        # First, generate all basic data from info
        permlist = info.deterministic_elemlist()
        cdef int order = len(permlist)
        permindices = {perm: i for i, perm in enumerate(permlist)} # 'invert' permlist
        identity = info.permgroup.identity
        # We use calloc to avoid calling free to unitialized pointers
        self.G = <qpms_finite_group_t *>calloc(1,sizeof(qpms_finite_group_t))
        if not self.G: raise MemoryError
        self.G[0].name = make_c_string(info.name)
        self.G[0].order = order
        self.G[0].idi = permindices[identity]
        self.G[0].mt = <qpms_gmi_t *>malloc(sizeof(qpms_gmi_t) * order * order)
        if not self.G[0].mt: raise MemoryError
        for i in range(order):
          for j in range(order):
            self.G[0].mt[i*order + j] = permindices[permlist[i] * permlist[j]]
        self.G[0].invi = <qpms_gmi_t *>malloc(sizeof(qpms_gmi_t) * order)
        if not self.G[0].invi: raise MemoryError
        for i in range(order):
            self.G[0].invi[i] = permindices[permlist[i]**-1]
        self.G[0].ngens = len(info.permgroupgens)
        self.G[0].gens = <qpms_gmi_t *>malloc(sizeof(qpms_gmi_t) * self.G[0].ngens)
        if not self.G[0].gens: raise MemoryError
        for i in range(self.G[0].ngens):
            self.G[0].gens[i] = permindices[info.permgroupgens[i]]
        self.G[0].permrep = <char **>calloc(order, sizeof(char *))
        if not self.G[0].permrep: raise MemoryError
        for i in range(order):
            self.G[0].permrep[i] = make_c_string(str(permlist[i]))
            if not self.G[0].permrep[i]: raise MemoryError
        self.G[0].permrep_nelem = info.permgroup.degree
        if info.rep3d is not None:
            self.G[0].rep3d = <qpms_irot3_t *>malloc(order * sizeof(qpms_irot3_t))
            for i in range(order):
                self.G[0].rep3d[i] = info.rep3d[permlist[i]].qd
        self.G[0].nirreps = len(info.irreps)
        self.G[0].irreps = <qpms_finite_group_irrep_t *>calloc(self.G[0].nirreps, sizeof(qpms_finite_group_irrep_t))
        if not self.G[0].irreps: raise MemoryError
        cdef int dim
        for iri, irname in enumerate(sorted(info.irreps.keys())):
            irrep = info.irreps[irname]
            is1d = isinstance(irrep[identity], (int, float, complex))
            dim = 1 if is1d else irrep[identity].shape[0]
            self.G[0].irreps[iri].dim = dim
            self.G[0].irreps[iri].name = <char *>make_c_string(irname)
            if not self.G[0].irreps[iri].name: raise MemoryError
            self.G[0].irreps[iri].m = <cdouble *>malloc(dim*dim*sizeof(cdouble)*order)
            if not self.G[0].irreps[iri].m: raise MemoryError
            if is1d:
                for i in range(order):
                    self.G[0].irreps[iri].m[i] = irrep[permlist[i]]
            else:
                for i in range(order):
                    for row in range(dim):
                        for col in range(dim):
                            self.G[0].irreps[iri].m[i*dim*dim + row*dim + col] = irrep[permlist[i]][row,col]
        self.G[0].elemlabels = <char **> 0 # Elem labels not yet implemented
        self.owns_data = True
        
    def __dealloc__(self):
        cdef qpms_gmi_t order
        if self.owns_data:
            if self.G:
                order = self.G[0].order
                free(self.G[0].name)
                free(self.G[0].mt)
                free(self.G[0].invi)
                free(self.G[0].gens)
                if self.G[0].permrep:
                    for i in range(order): free(self.G[0].permrep[i])
                free(self.G[0].permrep)
                if self.G[0].elemlabels: # this is not even contructed right now
                    for i in range(order): free(self.G[0].elemlabels[i])
                if self.G[0].irreps:
                    for iri in range(self.G[0].nirreps):
                        free(self.G[0].irreps[iri].name)
                        free(self.G[0].irreps[iri].m)
                free(self.G[0].irreps)
            free(self.G)
            self.G = <qpms_finite_group_t *>0
            self.owns_data = False

    cdef qpms_finite_group_t *rawpointer(self):
        return self.G

cdef class FinitePointGroupElement:
    '''TODO'''
    cdef readonly FinitePointGroup G
    cdef readonly qpms_gmi_t gmi
    def __cinit__(self, FinitePointGroup G, qpms_gmi_t gmi):
        self.G = G
        self.gmi = gmi

cdef class Particle:
    '''
    Wrapper over the qpms_particle_t structure.
    '''
    cdef qpms_particle_t p
    cdef readonly CTMatrix t # We hold the reference to the T-matrix to ensure correct reference counting

    def __cinit__(Particle self, pos, CTMatrix t):
        if(len(pos)>=2 and len(pos) < 4):
            self.p.pos.x = pos[0]
            self.p.pos.y = pos[1]
            self.p.pos.z = pos[2] if len(pos)==3 else 0
        else:
            raise ValueError("Position argument has to contain 3 or 2 cartesian coordinates")
        self.t = t
        self.p.tmatrix = self.t.rawpointer()

    cdef qpms_particle_t *rawpointer(Particle self):
        '''Pointer to the qpms_particle_p structure.
        '''
        return &(self.p)
    property rawpointer:
        def __get__(self):
            return <uintptr_t> &(self.p)

    cdef qpms_particle_t cval(Particle self):
        '''Provides a copy for assigning in cython code'''
        return self.p

    property x:
        def __get__(self):
            return self.p.pos.x
        def __set__(self,x):
            self.p.pos.x = x
    property y:
        def __get__(self):
            return self.p.pos.y
        def __set__(self,y):
            self.p.pos.y = y
    property z:
        def __get__(self):
            return self.p.pos.z
        def __set__(self,z):
            self.p.pos.z = z
    property pos:
        def __get__(self):
            return (self.p.pos.x, self.p.pos.y, self.p.pos.z)
        def __set__(self, pos):
            if(len(pos)>=2 and len(pos) < 4):
                self.p.pos.x = pos[0]
                self.p.pos.y = pos[1]
                self.p.pos.z = pos[2] if len(pos)==3 else 0
            else:
                raise ValueError("Position argument has to contain 3 or 2 cartesian coordinates")

cpdef void scatsystem_set_nthreads(long n):
    qpms_scatsystem_set_nthreads(n)
    return

cdef class ScatteringSystem:
    '''
    Wrapper over the C qpms_scatsys_t structure.
    '''
    cdef list basespecs # Here we keep the references to occuring basespecs
    #cdef list Tmatrices # Here we keep the references to occuring T-matrices
    cdef qpms_scatsys_t *s

    def __cinit__(self, particles, FinitePointGroup sym):
        '''TODO doc.
        Takes the particles (which have to be a sequence of instances of Particle),
        fills them together with their t-matrices to the "proto-qpms_scatsys_t"
        orig and calls qpms_scatsys_apply_symmetry
        (and then cleans orig)
        '''
        cdef qpms_scatsys_t orig # This should be automatically init'd to 0 (CHECKME)
        cdef qpms_ss_pi_t p_count = len(particles)
        cdef qpms_ss_tmi_t tm_count = 0
        tmindices = dict()
        tmobjs = list()
        self.basespecs=list()
        for p in particles: # find and enumerate unique t-matrices
            if id(p.t) not in tmindices:
                tmindices[id(p.t)] = tm_count
                tmobjs.append(p.t)
                tm_count += 1
        orig.tm_count = tm_count
        orig.p_count = p_count
        for tm in tmobjs: # create references to BaseSpec objects
            self.basespecs.append(tm.spec)
        try:
            orig.tm = <qpms_tmatrix_t **>malloc(orig.tm_count * sizeof(orig.tm[0]))
            if not orig.tm: raise MemoryError
            orig.p = <qpms_particle_tid_t *>malloc(orig.p_count * sizeof(orig.p[0]))
            if not orig.p: raise MemoryError
            for tmi in range(tm_count):
                orig.tm[tmi] = (<CTMatrix?>(tmobjs[tmi])).rawpointer()
            for pi in range(p_count):
                orig.p[pi].pos = (<Particle?>(particles[pi])).cval().pos
                orig.p[pi].tmatrix_id = tmindices[id(particles[pi].t)]
            self.s = qpms_scatsys_apply_symmetry(&orig, sym.rawpointer())
        finally:
            free(orig.tm)
            free(orig.p)

    def __dealloc__(self):
        qpms_scatsys_free(self.s)

    def particles_tmi(self):
        r = list()
        cdef qpms_ss_pi_t pi
        for pi in range(self.s[0].p_count):
            r.append(self.s[0].p[pi])
        return r

    property fecv_size: 
        def __get__(self): return self.s[0].fecv_size
    property saecv_sizes: 
        def __get__(self): 
            return [self.s[0].saecv_sizes[i] 
                for i in range(self.s[0].sym[0].nirreps)]
    property irrep_names: 
        def __get__(self): 
            return [string_c2py(self.s[0].sym[0].irreps[iri].name) 
                    if (self.s[0].sym[0].irreps[iri].name) else None
                for iri in range(self.s[0].sym[0].nirreps)]
    property nirreps: 
        def __get__(self): return self.s[0].sym[0].nirreps

    def pack_vector(self, vect, iri):
        if len(vect) != self.fecv_size: 
            raise ValueError("Length of a full vector has to be %d, not %d" 
                    % (self.fecv_size, len(vect)))
        vect = np.array(vect, dtype=complex, copy=False, order='C')
        cdef cdouble[::1] vect_view = vect;
        cdef np.ndarray[np.complex_t, ndim=1] target_np = np.empty(
                (self.saecv_sizes[iri],), dtype=complex, order='C')
        cdef cdouble[::1] target_view = target_np
        qpms_scatsys_irrep_pack_vector(&target_view[0], &vect_view[0], self.s, iri)
        return target_np
    def unpack_vector(self, packed, iri):
        if len(packed) != self.saecv_sizes[iri]: 
            raise ValueError("Length of %d. irrep-packed vector has to be %d, not %d"
                    % (iri, self.saecv_sizes, len(packed)))
        packed = np.array(packed, dtype=complex, copy=False, order='C')
        cdef cdouble[::1] packed_view = packed
        cdef np.ndarray[np.complex_t, ndim=1] target_np = np.empty(
                (self.fecv_size,), dtype=complex)
        cdef cdouble[::1] target_view = target_np
        qpms_scatsys_irrep_unpack_vector(&target_view[0], &packed_view[0], 
                self.s, iri, 0)
        return target_np
    def pack_matrix(self, fullmatrix, iri):
        cdef size_t flen = self.s[0].fecv_size
        cdef size_t rlen = self.saecv_sizes[iri]
        fullmatrix = np.array(fullmatrix, dtype=complex, copy=False, order='C')
        if fullmatrix.shape != (flen, flen):
            raise ValueError("Full matrix shape should be (%d,%d), is %s."
                    % (flen, flen, repr(fullmatrix.shape)))
        cdef cdouble[:,::1] fullmatrix_view = fullmatrix
        cdef np.ndarray[np.complex_t, ndim=2] target_np = np.empty(
                (rlen, rlen), dtype=complex, order='C')
        cdef cdouble[:,::1] target_view = target_np
        qpms_scatsys_irrep_pack_matrix(&target_view[0][0], &fullmatrix_view[0][0],
                self.s, iri)
        return target_np
    def unpack_matrix(self, packedmatrix, iri):
        cdef size_t flen = self.s[0].fecv_size
        cdef size_t rlen = self.saecv_sizes[iri]
        packedmatrix = np.array(packedmatrix, dtype=complex, copy=False, order='C')
        if packedmatrix.shape != (rlen, rlen):
            raise ValueError("Packed matrix shape should be (%d,%d), is %s."
                    % (rlen, rlen, repr(packedmatrix.shape)))
        cdef cdouble[:,::1] packedmatrix_view = packedmatrix
        cdef np.ndarray[np.complex_t, ndim=2] target_np = np.empty(
                (flen, flen), dtype=complex, order='C')
        cdef cdouble[:,::1] target_view = target_np
        qpms_scatsys_irrep_unpack_matrix(&target_view[0][0], &packedmatrix_view[0][0],
                self.s, iri, 0)
        return target_np

    def modeproblem_matrix_full(self, double k):
        cdef size_t flen = self.s[0].fecv_size
        cdef np.ndarray[np.complex_t, ndim=2] target = np.empty(
                (flen,flen),dtype=complex, order='C')
        cdef cdouble[:,::1] target_view = target
        qpms_scatsys_build_modeproblem_matrix_full(&target_view[0][0], self.s, k)
        return target

    def modeproblem_matrix_packed(self, double k, qpms_iri_t iri, version='pR'):
        cdef size_t rlen = self.saecv_sizes[iri]
        cdef np.ndarray[np.complex_t, ndim=2] target = np.empty(
                (rlen,rlen),dtype=complex, order='C')
        cdef cdouble[:,::1] target_view = target
        if (version == 'R'):
            qpms_scatsys_build_modeproblem_matrix_irrep_packed_orbitorderR(&target_view[0][0], self.s, iri, k)
        elif (version == 'pR'):
          with nogil:
            qpms_scatsys_build_modeproblem_matrix_irrep_packed_parallelR(&target_view[0][0], self.s, iri, k)
        else:
            qpms_scatsys_build_modeproblem_matrix_irrep_packed(&target_view[0][0], self.s, iri, k)
        return target

    def translation_matrix_full(self, double k, J = QPMS_HANKEL_PLUS):
        cdef size_t flen = self.s[0].fecv_size
        cdef np.ndarray[np.complex_t, ndim=2] target = np.empty(
                (flen,flen),dtype=complex, order='C')
        cdef cdouble[:,::1] target_view = target
        qpms_scatsys_build_translation_matrix_e_full(&target_view[0][0], self.s, k, J)
        return target

    def translation_matrix_packed(self, double k, qpms_iri_t iri, J = QPMS_HANKEL_PLUS):
        cdef size_t rlen = self.saecv_sizes[iri]
        cdef np.ndarray[np.complex_t, ndim=2] target = np.empty(
                (rlen,rlen),dtype=complex, order='C')
        cdef cdouble[:,::1] target_view = target
        qpms_scatsys_build_translation_matrix_e_irrep_packed(&target_view[0][0],
                self.s, iri, k, J)
        return target
    
    def fullvec_psizes(self):
        cdef np.ndarray[int32_t, ndim=1] ar = np.empty((self.s[0].p_count,), dtype=np.int32)
        cdef int32_t[::1] ar_view = ar
        for pi in range(self.s[0].p_count):
            ar_view[pi] = self.s[0].tm[self.s[0].p[pi].tmatrix_id].spec[0].n
        return ar


    def fullvec_poffsets(self):
        cdef np.ndarray[intptr_t, ndim=1] ar = np.empty((self.s[0].p_count,), dtype=np.intp)
        cdef intptr_t[::1] ar_view = ar
        cdef intptr_t offset = 0
        for pi in range(self.s[0].p_count):
            ar_view[pi] = offset
            offset += self.s[0].tm[self.s[0].p[pi].tmatrix_id].spec[0].n
        return ar

    def positions(self):
        cdef np.ndarray[np.double_t, ndim=2] ar = np.empty((self.s[0].p_count, 3), dtype=float)
        cdef np.double_t[:,::1] ar_view = ar
        for pi in range(self.s[0].p_count):
            ar_view[pi,0] = self.s[0].p[pi].pos.x
            ar_view[pi,1] = self.s[0].p[pi].pos.y
            ar_view[pi,2] = self.s[0].p[pi].pos.z
        return ar
   
    def planewave_full(self, k_cart, E_cart):
        k_cart = np.array(k_cart)
        E_cart = np.array(E_cart)
        if k_cart.shape != (3,) or E_cart.shape != (3,):
            raise ValueError("k_cart and E_cart must be ndarrays of shape (3,)")
        cdef qpms_incfield_planewave_params_t p
        p.use_cartesian = 1
        p.k.cart.x = <cdouble>k_cart[0]
        p.k.cart.y = <cdouble>k_cart[1]
        p.k.cart.z = <cdouble>k_cart[2]
        p.E.cart.x = <cdouble>E_cart[0]
        p.E.cart.y = <cdouble>E_cart[1]
        p.E.cart.z = <cdouble>E_cart[2]
        cdef np.ndarray[np.complex_t, ndim=1] target_np = np.empty(
                (self.fecv_size,), dtype=complex)
        cdef cdouble[::1] target_view = target_np
        qpms_scatsys_incident_field_vector_full(&target_view[0],
                self.s, qpms_incfield_planewave, <void *>&p, 0)
        return target_np

    def apply_Tmatrices_full(self, a):
        if len(a) != self.fecv_size: 
            raise ValueError("Length of a full vector has to be %d, not %d" 
                    % (self.fecv_size, len(a)))
        a = np.array(a, dtype=complex, copy=False, order='C')
        cdef cdouble[::1] a_view = a;
        cdef np.ndarray[np.complex_t, ndim=1] target_np = np.empty(
                (self.fecv_size,), dtype=complex, order='C')
        cdef cdouble[::1] target_view = target_np
        qpms_scatsys_apply_Tmatrices_full(&target_view[0], &a_view[0], self.s)
        return target_np
    
    cdef qpms_scatsys_t *rawpointer(self):
        return self.s

    def scatter_solver(self, double k, iri=None):
        return ScatteringMatrix(self, k, iri)

cdef class ScatteringMatrix:
    '''
    Wrapper over the C qpms_ss_LU structure that keeps the factorised mode problem matrix.
    '''
    cdef ScatteringSystem ss # Here we keep the reference to the parent scattering system
    cdef qpms_ss_LU lu

    def __cinit__(self, ScatteringSystem ss, double k, iri=None):
        self.ss = ss
        # TODO? pre-allocate the matrix with numpy to make it transparent?
        if iri is None:
            self.lu = qpms_scatsys_build_modeproblem_matrix_full_LU(
                    NULL, NULL, ss.rawpointer(), k)
        else:
            self.lu = qpms_scatsys_build_modeproblem_matrix_irrep_packed_LU(
                    NULL, NULL, ss.rawpointer(), iri, k)

    def __dealloc__(self):
        qpms_ss_LU_free(self.lu)

    property iri:
        def __get__(self):
            return None if self.lu.full else self.lu.iri

    def __call__(self, a_inc):
        cdef size_t vlen
        cdef qpms_iri_t iri = -1;
        if self.lu.full:
            vlen = self.lu.ss[0].fecv_size
            if len(a_inc) != vlen:
                raise ValueError("Length of a full coefficient vector has to be %d, not %d"
                        % (vlen, len(a_inc)))
        else:
            iri = self.lu.iri
            vlen = self.lu.ss[0].saecv_sizes[iri]
            if len(a_inc) != vlen:
                raise ValueError("Length of a %d. irrep packed coefficient vector has to be %d, not %d"
                        % (iri, vlen, len(a_inc)))
        a_inc = np.array(a_inc, dtype=complex, copy=False, order='C')
        cdef const cdouble[::1] a_view = a_inc;
        cdef np.ndarray f = np.empty((vlen,), dtype=complex, order='C')
        cdef cdouble[::1] f_view = f
        qpms_scatsys_scatter_solve(&f_view[0], &a_view[0], self.lu)
        return f


