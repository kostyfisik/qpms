"""@package qpms_c
Cythonized parts of QPMS; mostly wrappers over the C data structures
to make them available in Python.
"""

# Cythonized parts of QPMS here
# -----------------------------

import numpy as np
from .qpms_cdefs cimport *
from .cyquaternions cimport IRot3, CQuat
from .cybspec cimport BaseSpec
from .cycommon cimport make_c_string
from .cycommon import string_c2py, PointGroupClass
from .cytmatrices cimport CTMatrix, TMatrixFunction, TMatrixGenerator
from .cymaterials cimport EpsMuGenerator
from libc.stdlib cimport malloc, free, calloc
import warnings


# Set custom GSL error handler. N.B. this is obviously not thread-safe.
cdef char *pgsl_err_reason
cdef char *pgsl_err_file
cdef int pgsl_err_line
cdef int pgsl_errno = 0
cdef int *pgsl_errno_ignorelist = NULL # list of ignored error codes, terminated by zero

# This error handler only sets the variables above
cdef void pgsl_error_handler(const char *reason, const char *_file, const int line, const int gsl_errno):
    global pgsl_err_reason, pgsl_err_file, pgsl_err_line, pgsl_errno, pgsl_errno_ignorelist
    cdef size_t i
    if(pgsl_errno_ignorelist):
        i = 0
        while pgsl_errno_ignorelist[i] != 0:
            if gsl_errno == pgsl_errno_ignorelist[i]:
                return
            i += 1
    pgsl_err_file = _file
    pgsl_err_reason = reason
    pgsl_errno = gsl_errno
    pgsl_err_line = line
    return

cdef const int* pgsl_set_ignorelist(const int *new_ignorelist):
    global pgsl_errno_ignorelist
    cdef const int *oldlist = pgsl_errno_ignorelist
    pgsl_errno_ignorelist = new_ignorelist
    return oldlist

cdef class pgsl_ignore_error():
    '''Context manager for setting a temporary list of errnos ignored by pgsl_error_handler.
    Always sets pgsl_error_handler.

    Performs pgsl_check_err() on exit unless 
    '''
    cdef const int *ignorelist_old
    cdef gsl_error_handler_t *old_handler
    cdef bint final_check
    cdef object ignorelist_python

    cdef int *ignorelist
    def __cinit__(self, *ignorelist, **kwargs):
        self.ignorelist = <int*>calloc((len(ignorelist)+1), sizeof(int))
        self.ignorelist_python = ignorelist
        for i in range(len(ignorelist)):
            self.ignorelist[i] = ignorelist[i]
        if "final_check" in kwargs.keys() and not kwargs["final_check"]:
            final_check = True
        final_check = False
    
    def __enter__(self):
        global pgsl_error_handler
        self.ignorelist_old = pgsl_set_ignorelist(self.ignorelist)
        self.old_handler = gsl_set_error_handler(pgsl_error_handler)
        return

    def __exit__(self, type, value, traceback):
        global pgsl_errno_ignorelist, pgsl_error_handler
        pgsl_set_ignorelist(self.ignorelist_old)
        gsl_set_error_handler(self.old_handler)
        if self.final_check:
            pgsl_check_err(retval = None, ignore = self.ignorelist_python)

    def __dealloc__(self):
        free(self.ignorelist)

def pgsl_check_err(retval = None, ignorelist = None):
    global pgsl_err_reason, pgsl_err_file, pgsl_err_line, pgsl_errno
    '''Check for possible errors encountered by pgsl_error_handler.
    Takes return value of a function as an optional argument, which is now ignored.
    '''
    cdef int errno_was
    if (pgsl_errno != 0):
        errno_was = pgsl_errno
        pgsl_errno = 0
        raise RuntimeError("Error %d in GSL calculation in %s:%d: %s" % (errno_was,
            string_c2py(pgsl_err_file), pgsl_err_line, string_c2py(pgsl_err_reason)))
    if (retval is not None and retval != 0 and ignorelist is not None and retval not in ignorelist):
        warnings.warn("Got non-zero return value %d" % retval)
    if retval is not None:
        return retval
    else:
        return 0

def set_gsl_pythonic_error_handling():
    '''
    Sets pgsl_error_handler as the GSL error handler to avoid crashing.
    '''
    gsl_set_error_handler(pgsl_error_handler)

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
    cdef readonly TMatrixFunction f # Reference to ensure correct reference counting


    def __cinit__(Particle self, pos, t, bspec = None):
        cdef TMatrixGenerator tgen
        cdef BaseSpec spec
        if(len(pos)>=2 and len(pos) < 4):
            self.p.pos.x = pos[0]
            self.p.pos.y = pos[1]
            self.p.pos.z = pos[2] if len(pos)==3 else 0
        else:
            raise ValueError("Position argument has to contain 3 or 2 cartesian coordinates")
        if isinstance(t, CTMatrix):
            tgen = TMatrixGenerator(t)
        elif isinstance(t, TMatrixGenerator):
            tgen = <TMatrixGenerator>t
        else: raise TypeError('t must be either CTMatrix or TMatrixGenerator, was %s' % str(type(t)))
        if bspec is not None:
            spec = bspec
        else:
            if isinstance(tgen.holder, CTMatrix):
                spec = (<CTMatrix>tgen.holder).spec
            else:
                raise ValueError("bspec argument must be specified separately for str(type(t))")
        self.f = TMatrixFunction(tgen, spec)
        self.p.tmg = self.f.rawpointer()
        # TODO non-trivial transformations later; if modified, do not forget to update ScatteringSystem constructor
        self.p.op = qpms_tmatrix_operation_noop

    def __dealloc__(self):
        qpms_tmatrix_operation_clear(&self.p.op)

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
    cdef list tmgobjs # here we keep the references to occuring TMatrixFunctions (and hence BaseSpecs and TMatrixGenerators)
    #cdef list Tmatrices # Here we keep the references to occuring T-matrices
    cdef EpsMuGenerator medium_holder # Here we keep the reference to medium generator
    cdef qpms_scatsys_t *s

    def check_s(self): # cdef instead?
        if self.s == <qpms_scatsys_t *>NULL:
            raise ValueError("ScatteringSystem's s-pointer not set. You must not use the default constructor; use the create() method instead")
        #TODO is there a way to disable the constructor outside this module?

    @staticmethod # We don't have any "standard" constructor for this right now
    def create(particles, medium, FinitePointGroup sym, cdouble omega): # TODO tolerances
        # These we are going to construct
        cdef ScatteringSystem self
        cdef _ScatteringSystemAtOmega pyssw
        
        cdef qpms_scatsys_t orig # This should be automatically init'd to 0 (CHECKME)
        cdef qpms_ss_pi_t pi, p_count = len(particles)
        cdef qpms_ss_tmi_t tmi, tm_count = 0
        cdef qpms_ss_tmgi_t tmgi, tmg_count = 0

        cdef qpms_scatsys_at_omega_t *ssw
        cdef qpms_scatsys_t *ss

        cdef Particle p
        
        tmgindices = dict()
        tmgobjs = list()
        tmindices = dict()
        tmlist = list()
        for p in particles: # find and enumerate unique t-matrix generators
            if p.p.op.typ != QPMS_TMATRIX_OPERATION_NOOP:
                raise NotImplementedError("currently, only no-op T-matrix operations are allowed in ScatteringSystem constructor")
            #tmg_key = id(p.f) # This causes a different generator for each particle -> SUPER SLOW
            tmg_key = (id(p.f.generator), id(p.f.spec))
            if tmg_key not in tmgindices:
                tmgindices[tmg_key] = tmg_count
                tmgobjs.append(p.f) # Save the references on BaseSpecs and TMatrixGenerators (via TMatrixFunctions)
                tmg_count += 1
            # Following lines have to be adjusted when nontrivial operations allowed:
            tm_derived_key = (tmg_key, None) # TODO unique representation of p.p.op instead of None
            if tm_derived_key not in tmindices:
                tmindices[tm_derived_key] = tm_count
                tmlist.append(tm_derived_key)
                tm_count += 1
        cdef EpsMuGenerator mediumgen = EpsMuGenerator(medium)
        orig.medium = mediumgen.g
        orig.tmg_count = tmg_count
        orig.tm_count = tm_count
        orig.p_count = p_count
        try:
            orig.tmg = <qpms_tmatrix_function_t *>malloc(orig.tmg_count * sizeof(orig.tmg[0]))
            if not orig.tmg: raise MemoryError
            orig.tm = <qpms_ss_derived_tmatrix_t *>malloc(orig.tm_count * sizeof(orig.tm[0]))
            if not orig.tm: raise MemoryError
            orig.p = <qpms_particle_tid_t *>malloc(orig.p_count * sizeof(orig.p[0]))
            if not orig.p: raise MemoryError
            for tmgi in range(orig.tmg_count):
                orig.tmg[tmgi] = (<TMatrixFunction?>tmgobjs[tmgi]).raw()
            for tmi in range(tm_count):
                tm_derived_key = tmlist[tmi]
                tmgi = tmgindices[tm_derived_key[0]]
                orig.tm[tmi].tmgi = tmgi
                orig.tm[tmi].op = qpms_tmatrix_operation_noop # TODO adjust when notrivial operations allowed
            for pi in range(p_count):
                p = particles[pi]
                tmg_key = (id(p.f.generator), id(p.f.spec))
                tm_derived_key = (tmg_key, None) # TODO unique representation of p.p.op instead of None
                orig.p[pi].pos = p.cval().pos
                orig.p[pi].tmatrix_id = tmindices[tm_derived_key]
            ssw = qpms_scatsys_apply_symmetry(&orig, sym.rawpointer(), omega, &QPMS_TOLERANCE_DEFAULT)
            ss = ssw[0].ss
        finally:
            free(orig.tmg)
            free(orig.tm)
            free(orig.p)
        self = ScatteringSystem()
        self.medium_holder = mediumgen
        self.s = ss
        self.tmgobjs = tmgobjs
        pyssw = _ScatteringSystemAtOmega()
        pyssw.ssw = ssw
        pyssw.ss_pyref = self
        return self, pyssw

    def __call__(self, cdouble omega):
        self.check_s()
        cdef _ScatteringSystemAtOmega pyssw = _ScatteringSystemAtOmega()
        pyssw.ssw = qpms_scatsys_at_omega(self.s, omega)
        pyssw.ss_pyref = self
        return pyssw

    def __dealloc__(self):
        if(self.s):
            qpms_scatsys_free(self.s)

    property particles_tmi:
      def __get__(self):
        self.check_s()
        r = list()
        cdef qpms_ss_pi_t pi
        for pi in range(self.s[0].p_count):
            r.append(self.s[0].p[pi])
        return r

    property fecv_size: 
        def __get__(self): 
            self.check_s()
            return self.s[0].fecv_size
    property saecv_sizes: 
        def __get__(self): 
            self.check_s()
            return [self.s[0].saecv_sizes[i] 
                for i in range(self.s[0].sym[0].nirreps)]
    property irrep_names: 
        def __get__(self): 
            self.check_s()
            return [string_c2py(self.s[0].sym[0].irreps[iri].name) 
                    if (self.s[0].sym[0].irreps[iri].name) else None
                for iri in range(self.s[0].sym[0].nirreps)]
    property nirreps: 
        def __get__(self): 
            self.check_s()
            return self.s[0].sym[0].nirreps

    def pack_vector(self, vect, iri):
        self.check_s()
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
        self.check_s()
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
        self.check_s()
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
        self.check_s()
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

    def translation_matrix_full(self, double k, J = QPMS_HANKEL_PLUS):
        self.check_s()
        cdef size_t flen = self.s[0].fecv_size
        cdef np.ndarray[np.complex_t, ndim=2] target = np.empty(
                (flen,flen),dtype=complex, order='C')
        cdef cdouble[:,::1] target_view = target
        qpms_scatsys_build_translation_matrix_e_full(&target_view[0][0], self.s, k, J)
        return target

    def translation_matrix_packed(self, double k, qpms_iri_t iri, J = QPMS_HANKEL_PLUS):
        self.check_s()
        cdef size_t rlen = self.saecv_sizes[iri]
        cdef np.ndarray[np.complex_t, ndim=2] target = np.empty(
                (rlen,rlen),dtype=complex, order='C')
        cdef cdouble[:,::1] target_view = target
        qpms_scatsys_build_translation_matrix_e_irrep_packed(&target_view[0][0],
                self.s, iri, k, J)
        return target
    
    property fullvec_psizes:
      def __get__(self):
        self.check_s()
        cdef np.ndarray[int32_t, ndim=1] ar = np.empty((self.s[0].p_count,), dtype=np.int32)
        cdef int32_t[::1] ar_view = ar
        for pi in range(self.s[0].p_count):
            ar_view[pi] = self.s[0].tm[self.s[0].p[pi].tmatrix_id].spec[0].n
        return ar


    property fullvec_poffsets:
      def __get__(self):
        self.check_s()
        cdef np.ndarray[intptr_t, ndim=1] ar = np.empty((self.s[0].p_count,), dtype=np.intp)
        cdef intptr_t[::1] ar_view = ar
        cdef intptr_t offset = 0
        for pi in range(self.s[0].p_count):
            ar_view[pi] = offset
            offset += self.s[0].tm[self.s[0].p[pi].tmatrix_id].spec[0].n
        return ar

    property positions:
      def __get__(self):
        self.check_s()
        cdef np.ndarray[np.double_t, ndim=2] ar = np.empty((self.s[0].p_count, 3), dtype=float)
        cdef np.double_t[:,::1] ar_view = ar
        for pi in range(self.s[0].p_count):
            ar_view[pi,0] = self.s[0].p[pi].pos.x
            ar_view[pi,1] = self.s[0].p[pi].pos.y
            ar_view[pi,2] = self.s[0].p[pi].pos.z
        return ar
   
    def planewave_full(self, k_cart, E_cart):
        self.check_s()
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

cdef class _ScatteringSystemAtOmega:
    '''
    Wrapper over the C qpms_scatsys_at_omega_t structure
    that keeps the T-matrix and background data evaluated
    at specific frequency.
    '''
    cdef qpms_scatsys_at_omega_t *ssw
    cdef ScatteringSystem ss_pyref
    
    def check(self): # cdef instead?
        if not self.ssw:
            raise ValueError("_ScatteringSystemAtOmega's ssw-pointer not set. You must not use the default constructor; ScatteringSystem.create() instead")
        self.ss_pyref.check_s()
        #TODO is there a way to disable the constructor outside this module?

    def __dealloc__(self):
        if (self.ssw):
            qpms_scatsys_at_omega_free(self.ssw)

    def apply_Tmatrices_full(self, a):
        self.check()
        if len(a) != self.fecv_size: 
            raise ValueError("Length of a full vector has to be %d, not %d" 
                    % (self.fecv_size, len(a)))
        a = np.array(a, dtype=complex, copy=False, order='C')
        cdef cdouble[::1] a_view = a;
        cdef np.ndarray[np.complex_t, ndim=1] target_np = np.empty(
                (self.fecv_size,), dtype=complex, order='C')
        cdef cdouble[::1] target_view = target_np
        qpms_scatsysw_apply_Tmatrices_full(&target_view[0], &a_view[0], self.ssw)
        return target_np
    
    cdef qpms_scatsys_at_omega_t *rawpointer(self):
        return self.ssw

    def scatter_solver(self, iri=None):
        self.check()
        return ScatteringMatrix(self, iri)

    property fecv_size: 
        def __get__(self): return self.ss_pyref.fecv_size
    property saecv_sizes: 
        def __get__(self): return self.ss_pyref.saecv_sizes
    property irrep_names: 
        def __get__(self): return self.ss_pyref.irrep_names
    property nirreps: 
        def __get__(self): return self.ss_pyref.nirreps

    def modeproblem_matrix_full(self):
        self.check()
        cdef size_t flen = self.ss_pyref.s[0].fecv_size
        cdef np.ndarray[np.complex_t, ndim=2] target = np.empty(
                (flen,flen),dtype=complex, order='C')
        cdef cdouble[:,::1] target_view = target
        qpms_scatsysw_build_modeproblem_matrix_full(&target_view[0][0], self.ssw)
        return target

    def modeproblem_matrix_packed(self, qpms_iri_t iri, version='pR'):
        self.check()
        cdef size_t rlen = self.saecv_sizes[iri]
        cdef np.ndarray[np.complex_t, ndim=2] target = np.empty(
                (rlen,rlen),dtype=complex, order='C')
        cdef cdouble[:,::1] target_view = target
        if (version == 'R'):
            qpms_scatsysw_build_modeproblem_matrix_irrep_packed_orbitorderR(&target_view[0][0], self.ssw, iri)
        elif (version == 'pR'):
          with nogil:
            qpms_scatsysw_build_modeproblem_matrix_irrep_packed(&target_view[0][0], self.ssw, iri)
        else:
            qpms_scatsysw_build_modeproblem_matrix_irrep_packed_serial(&target_view[0][0], self.ssw, iri)
        return target


cdef class ScatteringMatrix:
    '''
    Wrapper over the C qpms_ss_LU structure that keeps the factorised mode problem matrix.
    '''
    cdef _ScatteringSystemAtOmega ssw # Here we keep the reference to the parent scattering system
    cdef qpms_ss_LU lu

    def __cinit__(self, _ScatteringSystemAtOmega ssw, iri=None):
        ssw.check()
        self.ssw = ssw
        # TODO? pre-allocate the matrix with numpy to make it transparent?
        if iri is None:
            self.lu = qpms_scatsysw_build_modeproblem_matrix_full_LU(
                    NULL, NULL, ssw.rawpointer())
        else:
            self.lu = qpms_scatsysw_build_modeproblem_matrix_irrep_packed_LU(
                    NULL, NULL, ssw.rawpointer(), iri)

    def __dealloc__(self):
        qpms_ss_LU_free(self.lu)

    property iri:
        def __get__(self):
            return None if self.lu.full else self.lu.iri

    def __call__(self, a_inc):
        cdef size_t vlen
        cdef qpms_iri_t iri = -1;
        if self.lu.full:
            vlen = self.lu.ssw[0].ss[0].fecv_size
            if len(a_inc) != vlen:
                raise ValueError("Length of a full coefficient vector has to be %d, not %d"
                        % (vlen, len(a_inc)))
        else:
            iri = self.lu.iri
            vlen = self.lu.ssw[0].ss[0].saecv_sizes[iri]
            if len(a_inc) != vlen:
                raise ValueError("Length of a %d. irrep packed coefficient vector has to be %d, not %d"
                        % (iri, vlen, len(a_inc)))
        a_inc = np.array(a_inc, dtype=complex, copy=False, order='C')
        cdef const cdouble[::1] a_view = a_inc;
        cdef np.ndarray f = np.empty((vlen,), dtype=complex, order='C')
        cdef cdouble[::1] f_view = f
        qpms_scatsys_scatter_solve(&f_view[0], &a_view[0], self.lu)
        return f

def pitau(double theta, qpms_l_t lMax, double csphase = -1):
    if(abs(csphase) != 1):
        raise ValueError("csphase must be 1 or -1, is %g" % csphase)
    cdef size_t nelem = qpms_lMax2nelem(lMax)
    cdef np.ndarray[np.float_t, ndim=1] lega = np.empty((nelem,), dtype=float)
    cdef np.ndarray[np.float_t, ndim=1] pia = np.empty((nelem,), dtype=float)
    cdef np.ndarray[np.float_t, ndim=1] taua = np.empty((nelem,), dtype=float)
    cdef double[::1] leg = lega
    cdef double[::1] pi = pia
    cdef double[::1] tau = taua
    qpms_pitau_fill(&leg[0], &pi[0], &tau[0], theta, lMax, csphase)
    return (lega, pia, taua)

def linton_gamma(cdouble x):
    return clilgamma(x)

def linton_gamma_real(double x):
    return lilgamma(x)

def gamma_inc(double a, cdouble x, int m = 0):
    cdef qpms_csf_result res
    with pgsl_ignore_error(15): #15 is underflow
        complex_gamma_inc_e(a, x, m, &res)
    return (res.val, res.err)

def gamma_inc_series(double a, cdouble x):
    cdef qpms_csf_result res
    with pgsl_ignore_error(15): #15 is underflow
        cx_gamma_inc_series_e(a, x, &res)
    return (res.val, res.err)

def gamma_inc_CF(double a, cdouble x):
    cdef qpms_csf_result res
    with pgsl_ignore_error(15): #15 is underflow
        cx_gamma_inc_CF_e(a, x, &res)
    return (res.val, res.err)

def lll_reduce(basis, double delta=0.75):
    """
    Lattice basis reduction with the Lenstra-Lenstra-Lovász algorithm.

    basis is array_like with dimensions (n, d), where
    n is the size of the basis (dimensionality of the lattice)
    and d is the dimensionality of the space into which the lattice
    is embedded.
    """
    basis = np.array(basis, copy=True, order='C', dtype=np.double)
    if len(basis.shape) != 2:
        raise ValueError("Expected two-dimensional array (got %d-dimensional)"
                % len(basis.shape))
    cdef size_t n, d
    n, d = basis.shape
    if n > d:
        raise ValueError("Real space dimensionality (%d) cannot be smaller than"
                "the dimensionality of the lattice (%d) embedded into it."
                % (d, n))
    cdef double [:,:] basis_view = basis
    if 0 != qpms_reduce_lattice_basis(&basis_view[0,0], n, d, delta):
        raise RuntimeError("Something weird happened")
    return basis


cdef PGen get_PGen_direct(direct_basis, bint include_origin=False, double layers=30):
    dba = np.array(direct_basis)
    if not (dba.shape == (2,2)):
        raise NotImplementedError
    cdef cart2_t b1, b2
    b1.x = dba[0,0]
    b1.y = dba[0,1]
    b2.x = dba[1,0]
    b2.y = dba[0,1]
    cdef double maxR = layers*max(cart2norm(b1), cart2norm(b2))
    return PGen_xyWeb_new(b1, b2, BASIS_RTOL, CART2_ZERO, 0, include_origin, maxR, False)

cdef double get_unitcell_volume(direct_basis):
    dba = np.array(direct_basis)
    if not (dba.shape == (2,2)):
        raise NotImplementedError
    cdef cart2_t b1, b2
    b1.x = dba[0,0]
    b1.y = dba[0,1]
    b2.x = dba[1,0]
    b2.y = dba[0,1]
    return l2d_unitcell_area(b1, b2)

cdef PGen get_PGen_reciprocal2pi(direct_basis, double layers = 30):
    dba = np.array(direct_basis)
    if not (dba.shape == (2,2)):
        raise NotImplementedError
    cdef cart2_t b1, b2, rb1, rb2
    b1.x = dba[0,0]
    b1.y = dba[0,1]
    b2.x = dba[1,0]
    b2.y = dba[0,1]
    if(l2d_reciprocalBasis2pi(b1, b2, &rb1, &rb2) != 0):
            raise RuntimeError
    cdef double maxK = layers*max(cart2norm(rb1), cart2norm(rb2))
    return PGen_xyWeb_new(rb1, rb2, BASIS_RTOL, CART2_ZERO,
            0, True, maxK, False)

cdef class Ewald3Calculator:
    '''Wrapper class over qpms_ewald3_constants_t.
    
    Mainly for testing low-level scalar Ewald summation functionality.'''
    cdef qpms_ewald3_constants_t *c

    def __cinit__(self, qpms_l_t lMax, int csphase = -1):
        if (csphase != -1 and csphase != 1):
            raise ValueError("csphase must be +1 or -1, not %d" % csphase)
        self.c = qpms_ewald3_constants_init(lMax, csphase)

    def __dealloc__(self):
        qpms_ewald3_constants_free(self.c)

    def sigma0(self, double eta, cdouble wavenumber, do_err = False):
        cdef int retval
        cdef double err
        cdef cdouble result
        retval = ewald3_sigma0(&result, &err, self.c, eta, wavenumber)
        if retval:
            raise RuntimeError("ewald3_sigma0 returned non-zero value (%d)" % retval)
        if do_err:
            return (result, err)
        else:
            return result

    def sigma_short(self, double eta, cdouble wavenumber, direct_basis, wavevector, particle_shift, do_err=False):
        # FIXME now only 2d XY lattice in 3D is implemented here, we don't even do proper dimensionality checks.
        cdef cart3_t beta, pshift
        beta.x = wavevector[0]
        beta.y = wavevector[1]
        beta.z = 0
        pshift.x = particle_shift[0]
        pshift.y = particle_shift[1]
        pshift.z = 0
        cdef qpms_l_t n = self.c[0].nelem_sc
        cdef np.ndarray[complex, ndim=1] result = np.empty((n,), dtype=complex)
        cdef cdouble[::1] result_v = result
        cdef np.ndarray[double, ndim=1] err
        cdef double[::1] err_v
        if do_err:
            err = np.empty((n,), dtype=np.double)
            err_v = err
        cdef bint include_origin = not (particle_shift[0] == 0 and particle_shift[1] == 0)
        cdef PGen rgen = get_PGen_direct(direct_basis, include_origin)
        cdef int retval = ewald3_sigma_short(&result_v[0], &err_v[0] if do_err else NULL, 
                self.c, eta, wavenumber, LAT_2D_IN_3D_XYONLY, &rgen, False, beta, pshift)
        if rgen.stateData: PGen_destroy(&rgen)
        if retval: raise RuntimeError("ewald3_sigma_short returned %d" % retval)
        if do_err:
            return (result, err)
        else:
            return result

    def sigma_long(self, double eta, cdouble wavenumber, direct_basis, wavevector, particle_shift, do_err=False):
        # FIXME now only 2d XY lattice in 3D is implemented here, we don't even do proper dimensionality checks.
        cdef cart3_t beta, pshift
        beta.x = wavevector[0]
        beta.y = wavevector[1]
        beta.z = 0
        pshift.x = particle_shift[0]
        pshift.y = particle_shift[1]
        pshift.z = 0
        cdef qpms_l_t n = self.c[0].nelem_sc
        cdef np.ndarray[complex, ndim=1] result = np.empty((n,), dtype=complex)
        cdef cdouble[::1] result_v = result
        cdef np.ndarray[double, ndim=1] err
        cdef double[::1] err_v
        if do_err:
            err = np.empty((n,), dtype=np.double)
            err_v = err
        cdef PGen kgen = get_PGen_reciprocal2pi(direct_basis)
        cdef double unitcell_volume = get_unitcell_volume(direct_basis)
        cdef int retval = ewald3_sigma_long(&result_v[0], &err_v[0] if do_err else NULL,
                self.c, eta, wavenumber, unitcell_volume, LAT_2D_IN_3D_XYONLY, &kgen, False, beta, pshift)

        if kgen.stateData: PGen_destroy(&kgen)
        if retval: raise RuntimeError("ewald3_sigma_long returned %d" % retval)
        if do_err:
            return (result, err)
        else:
            return result
        



