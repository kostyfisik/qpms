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
from .cycommon import string_c2py, PointGroupClass, BesselType
from .cytmatrices cimport CTMatrix, TMatrixFunction, TMatrixGenerator, TMatrixInterpolator
from .cymaterials cimport EpsMuGenerator, EpsMu
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

    property order: # LPTODO might instead be __len__() if iterable at some point
        def __get__(self):
            return self.G[0].order

    @staticmethod 
    def TRIVIAL(): # quite ugly
        from .symmetries import point_group_info
        return FinitePointGroup(point_group_info['trivial_g'])


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
        """Particle object constructor.

        Parameters
        ----------
        pos : (x, y, z) 
            Particle position in cartesian coordinates

        t : TMatrixGenerator or CTMatrix or TMatrixInterpolator
            T-matrix specification for the particle.

        bspec : BaseSpec
            WSWF basis specification for the particle. Might be omitted if 
            t is a CTMatrix.
        """
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
        elif isinstance(t, TMatrixInterpolator):
            tgen = TMatrixGenerator(t)
            warnings.warn("Initialising a particle with interpolated T-matrix values. Imaginary frequencies will be discarded and mode search algorithm will yield nonsense (just saying).")
        elif isinstance(t, TMatrixGenerator):
            tgen = <TMatrixGenerator>t
        else: raise TypeError('t must be either CTMatrix or TMatrixGenerator, was %s' % str(type(t)))
        if bspec is not None:
            spec = bspec
        else:
            if isinstance(tgen.holder, CTMatrix):
                spec = (<CTMatrix>tgen.holder).spec
            else:
                raise ValueError("bspec argument must be specified separately for %s" % str(type(t)))
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
    Wrapper over the C qpms_scatsys_t structure, representing a collection
    of scatterers, finite or periodic.

    Currently, it does not have a standard constructor. Use the
    ScatteringSystem.create() method instead.
    '''
    cdef list tmgobjs # here we keep the references to occuring TMatrixFunctions (and hence BaseSpecs and TMatrixGenerators)
    #cdef list Tmatrices # Here we keep the references to occuring T-matrices
    cdef EpsMuGenerator medium_holder # Here we keep the reference to medium generator
    cdef qpms_scatsys_t *s
    cdef FinitePointGroup sym

    cdef qpms_iri_t iri_py2c(self, iri, allow_None = True):
        if iri is None and allow_None:
            return QPMS_NO_IRREP
        cdef qpms_iri_t nir = self.nirreps
        cdef qpms_iri_t ciri = iri
        if ciri < 0 or ciri > nir:
            raise ValueError("Invalid irrep index %s (of %d irreps)", str(iri), self.nirreps)
        return ciri

    def check_s(self): # cdef instead?
        if self.s == <qpms_scatsys_t *>NULL:
            raise ValueError("ScatteringSystem's s-pointer not set. You must not use the default constructor; use the create() method instead")
        #TODO is there a way to disable the constructor outside this module?

    @staticmethod # We don't have any "standard" constructor for this right now
    def create(particles, medium, cdouble omega, FinitePointGroup sym = FinitePointGroup.TRIVIAL(),
            latticebasis = None): # TODO tolerances
        """(Non-standard) constructor of ScatteringSystem

        Parameters
        ----------
        particles : list of Particles objects
            Scatterers to be included in the system. These scatterers are then
            copied around by the constructor using the point group symmetry defined in sym.
        medium : EpsMu or EpsMuGenerator
            Material properties of the background medium.
        omega : complex
            Any valid angular frequency for the initial T-matrix evaluation.
            It must be a value which all the T-matrix generators and interpolators 
            referenced by particles can evaluate.
        sym : FinitePointGroup
            Symmetry group for a finite system.
            Defaults to the trivial point group FinitePointGroup.TRIVIAL().
            Currently, this must be left trivial for periodic systems.
        latticebasis : None or array_like of float type.
            Lattice base vectors of a periodic system in cartesian coordinates.
            The shape must be [d][3], where d = 1, 2 or 3 determines the lattice dimension.
            Defaults to None, i.e. a finite system.

        Returns
        -------
        ss : ScatteringSystem
            Object representing the system of compact scatterers.
        ssw : _ScatteringSystemAtOmega
            Object representing ss evaluated at omega (angular frequency used for the 
            initial evaluation).

        Note
        ----
        Currently, the ScatterinSystem's T-matrices need to be evaluated for at least
        one frequency in order to initialise the symmetry structure. For this reason,
        this "constructor" creates an instance of ScatteringSystem and an instance
        of _ScatteringSystemAtOmega at the same time, so that the evaluated T-matrices
        can be used right away if needed. This is why this approach is used instead
        of the usual __init__ method.
        """

        if latticebasis is not None and sym.order != 1:
            raise NotImplementedError("Periodic systems don't currently support nontrivial point group symmetries")

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
            if latticebasis is not None: # periodic system
                assert(len(latticebasis) <= 3 and len(latticebasis) > 0)
                orig.lattice_dimension = len(latticebasis)
                for d in range(len(latticebasis)):
                    orig.per.lattice_basis[d] = {'x' : latticebasis[d][0], 'y' : latticebasis[d][1], 'z' : latticebasis[d][2] if len(latticebasis[d]) >= 3 else 0}
            else: orig.lattice_dimension = 0
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
        self.sym = sym
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
        """Length of the full excitation coefficient vector"""
        def __get__(self): 
            self.check_s()
            return self.s[0].fecv_size
    property saecv_sizes:
        """Lengths of the partial symmetry-adapted excitation coefficient vectors"""
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
        """Number of irreducible representations of the scattering system point symmetry group"""
        def __get__(self): 
            self.check_s()
            return self.s[0].sym[0].nirreps
    property lattice_dimension:
        def __get__(self):
            return self.s[0].lattice_dimension

    property unitcell_volume:
        def __get__(self):
            self.check_s()
            if self.lattice_dimension:
                return self.s[0].per.unitcell_volume
            else:
                return None

    property eta:
        """Ewald parameter η"""
        def __get__(self):
            self.check_s()
            if self.lattice_dimension:
                return self.s[0].per.eta
            else:
                return None

        def __set__(self, eta):
            self.check_s()
            if self.lattice_dimension:
                self.s[0].per.eta = eta
            else:
                raise AttributeError("Cannot set Ewald parameter for finite system") # different exception?


    def pack_vector(self, vect, iri):
        """Converts (projects) a full excitation coefficient vector into an irrep subspace.

        Parameters
        ----------
        vect : array_like of shape (self.fecv_size,)
            The full excitation coefficient vector to be converted.
        iri : int
            Index of the irreducible representation.

        Returns
        -------
        packed : ndarray of shape (self.saecv_sizes[iri],)
            "irrep-packed" excitation vector: the part of vect belonging to the 
            irrep indexed by iri, in terms of given irrep's internal coordinates.

        See Also
        --------
        unpack_vector : The beckward (not exactly inverse) operation.
        pack_matrix : Corresponding conversion for matrices.
        """
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
        """Unpacks an "irrep-packed" excitation coefficient vector to full coordinates.

        Parameters
        ----------
        packed : array_like of shape (self.saecv_sizes[iri],)
            Excitation coefficient vector part belonging to the irreducible representation
            indexed by iri, in terms of given irrep's internal coordinates.
        iri : int
            Index of the irreducible representation.

        Returns
        -------
        vect : ndarray of shape (self.fecv_size,)
            The contribution to a full excitation coefficient vector from iri'th irrep.

        See Also
        --------
        pack_vector : The inverse operation.
        unpack_matrix : Corresponding conversion for matrices.
        """
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
        """Converts (projects) a matrix into an irrep subspace.

        Parameters
        ----------
        fullmatrix : array_like of shape (self.fecv_size, self.fecv_size)
            The full matrix (operating on the exctitation coefficient vectors) to be converted.
        iri : int
            Index of the irreducible representation.

        Returns
        -------
        packedmatrix : ndarray of shape (self.saecv_sizes[iri], self.saecv_sizes[iri])
            "irrep-packed" matrix: the part of fullmatrix belonging to the 
            irrep indexed by iri, in terms of given irrep's internal coordinates.

        See Also
        --------
        unpack_matrix : The beckward (not exactly inverse) operation.
        pack_vector : Corresponding conversion for excitation coefficient vectors.
        """
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
        """Unpacks an "irrep-packed" excitation coefficient vector to full coordinates.

        Parameters
        ----------
        packedmatrix : array_like of shape (self.saecv_sizes[iri], self.saecv_sizes[iri])
            A matrix (operating on the excitation coefficient vectors) part belonging to the 
            irreducible representation indexed by iri, in terms of given irrep's internal
            coordinates.
        iri : int
            Index of the irreducible representation.

        Returns
        -------
        fullmatrix : ndarray of shape (self.fecv_size, self.fecv_size)
            The iri'th irrep contribution to a full matrix. 

        See Also
        --------
        pack_matrix : The inverse operation.
        unpack_vector : Corresponding conversion for excitation coefficient vectors.
        """
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

    def translation_matrix_full(self, cdouble wavenumber, blochvector = None, J = QPMS_HANKEL_PLUS):
        """Constructs full translation matrix of a scattering system.

        This method enables to use any wave number for the background medium ignoring the
        background EpsMuGenerator), using only system's particle positions (and lattice
        basis for infinite system).

        Parameters
        ----------
        wavenumber : complex
            Wave number of the medium

        blochvector : array_like of shape (3,)
            Bloch vector (only for periodic systems)

        J : BesselType
            Optionally, one can replace Hankel functions of the first kind with different
            Bessel functions.
        
        See Also
        --------
        ScatteringSystemAtOmega.translation_matrix_full : Translation matrix at a given frequency.
        """
        self.check_s()
        cdef size_t flen = self.s[0].fecv_size
        cdef np.ndarray[np.complex_t, ndim=2] target = np.empty(
                (flen,flen),dtype=complex, order='C')
        cdef cdouble[:,::1] target_view = target
        cdef cart3_t blochvector_c
        if self.lattice_dimension == 0:
            if blochvector is None:
                qpms_scatsys_build_translation_matrix_e_full(&target_view[0][0], self.s, wavenumber, J)
            else: raise ValueError("Can't use blochvector with non-periodic system")
        else:
            if blochvector is None: raise ValueError("Valid blochvector must be specified for periodic system")
            else:
                if J != QPMS_HANKEL_PLUS:
                    raise NotImplementedError("Translation operators based on other than Hankel+ functions not supperted in periodic systems")
                blochvector_c = {'x': blochvector[0], 'y': blochvector[1], 'z': blochvector[2]}
                with pgsl_ignore_error(15):
                    qpms_scatsys_periodic_build_translation_matrix_full(&target_view[0][0], self.s, wavenumber, &blochvector_c)
        return target

    def translation_matrix_packed(self, cdouble wavenumber, qpms_iri_t iri, J = QPMS_HANKEL_PLUS):
        self.check_s()
        cdef size_t rlen = self.saecv_sizes[iri]
        cdef np.ndarray[np.complex_t, ndim=2] target = np.empty(
                (rlen,rlen),dtype=complex, order='C')
        cdef cdouble[:,::1] target_view = target
        qpms_scatsys_build_translation_matrix_e_irrep_packed(&target_view[0][0],
                self.s, iri, wavenumber, J)
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
            offset += qpms_ss_bspec_pi(self.s, pi)[0].n
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

    def find_modes(self, cdouble omega_centre, double omega_rr, double omega_ri, iri = None,
            blochvector = None,
            size_t contour_points = 20, double rank_tol = 1e-4, size_t rank_min_sel=1,
            double res_tol = 0):
        """
        Attempts to find the eigenvalues and eigenvectors using Beyn's algorithm.
        """

        cdef beyn_result_t *res
        cdef double blochvector_c[3]
        if self.lattice_dimension == 0:
            assert(blochvector is None)
            res = qpms_scatsys_finite_find_eigenmodes(self.s, self.iri_py2c(iri),
                omega_centre, omega_rr, omega_ri, contour_points, 
                rank_tol, rank_min_sel, res_tol)
        else:
            if iri is not None: raise NotImplementedError("Irreps decomposition not yet supported for periodic systems")
            blochvector_c[0] = blochvector[0]
            blochvector_c[1] = blochvector[1]
            blochvector_c[2] = blochvector[2] if len(blochvector) > 2 else 0
            res = qpms_scatsys_periodic_find_eigenmodes(self.s, blochvector_c,
                    omega_centre, omega_rr, omega_ri, contour_points, rank_tol, rank_min_sel, res_tol)
        if res == NULL: raise RuntimeError

        cdef size_t neig = res[0].neig, i, j
        cdef size_t vlen = res[0].vlen # should be equal to self.s.fecv_size

        cdef np.ndarray[complex, ndim=1] eigval = np.empty((neig,), dtype=complex)
        cdef cdouble[::1] eigval_v = eigval
        cdef np.ndarray[complex, ndim=1] eigval_err = np.empty((neig,), dtype=complex)
        cdef cdouble[::1] eigval_err_v = eigval_err
        cdef np.ndarray[double, ndim=1] residuals = np.empty((neig,), dtype=np.double)
        cdef double[::1] residuals_v = residuals
        cdef np.ndarray[complex, ndim=2] eigvec = np.empty((neig,vlen),dtype=complex)
        cdef cdouble[:,::1] eigvec_v = eigvec
        cdef np.ndarray[double, ndim=1] ranktest_SV = np.empty((vlen), dtype=np.double)
        cdef double[::1] ranktest_SV_v = ranktest_SV

        for i in range(neig):
            eigval_v[i] = res[0].eigval[i]
            eigval_err_v[i] = res[0].eigval_err[i]
            residuals_v[i] = res[0].residuals[i]
            for j in range(vlen):
                eigvec_v[i,j] = res[0].eigvec[i*vlen + j]
        for i in range(vlen):
            ranktest_SV_v[i] = res[0].ranktest_SV[i]

        zdist = eigval - omega_centre
        eigval_inside_metric = np.hypot(zdist.real / omega_rr, zdist.imag / omega_ri)

        beyn_result_free(res)
        retdict = {
                'eigval':eigval,
                'eigval_inside_metric':eigval_inside_metric,
                'eigvec':eigvec,
                'residuals':residuals,
                'eigval_err':eigval_err,
                'ranktest_SV':ranktest_SV,
                'iri': iri,
                'blochvector': blochvector
        }

        return retdict

    def scattered_E(self, cdouble wavenumber, scatcoeffvector_full, evalpos, bint alt=False, btyp=BesselType.HANKEL_PLUS):
        cdef qpms_bessel_t btyp_c = BesselType(btyp)
        evalpos = np.array(evalpos, dtype=float, copy=False)
        if evalpos.shape[-1] != 3:
            raise ValueError("Last dimension of evalpos has to be 3")
        cdef np.ndarray[double,ndim=2] evalpos_a = evalpos.reshape(-1,3)
        cdef np.ndarray[dtype=complex, ndim=1] scv = np.array(scatcoeffvector_full, copy=False)
        cdef cdouble[::1] scv_view = scv
        cdef np.ndarray[complex, ndim=2] results = np.empty((evalpos_a.shape[0],3), dtype=complex)
        cdef ccart3_t res
        cdef cart3_t pos
        cdef size_t i
        for i in range(evalpos_a.shape[0]):
            pos.x = evalpos_a[i,0]
            pos.y = evalpos_a[i,1]
            pos.z = evalpos_a[i,2]
            if alt:
                res = qpms_scatsys_scattered_E__alt(self.s, btyp_c, wavenumber, &scv_view[0], pos)
            else:
                res = qpms_scatsys_scattered_E(self.s, btyp_c, wavenumber, &scv_view[0], pos)
            results[i,0] = res.x
            results[i,1] = res.y
            results[i,2] = res.z
        return results.reshape(evalpos.shape)

def empty_lattice_modes_xy(EpsMu epsmu, reciprocal_basis, wavevector, double maxomega):
    '''Empty (2D, xy-plane) lattice mode (diffraction order) frequencies of a non-dispersive medium.

    reciprocal_basis is of the "mutliplied by 2п" type.
    '''
    cdef double *omegas_c
    cdef size_t n
    cdef cart2_t k_, b1, b2
    k_.x = wavevector[0]
    k_.y = wavevector[1]
    b1.x = reciprocal_basis[0][0]
    b1.y = reciprocal_basis[0][1]
    b2.x = reciprocal_basis[1][0]
    b2.y = reciprocal_basis[1][1]
    if(epsmu.n.imag != 0):
        warnings.warn("Got complex refractive index", epsmu.n, "ignoring the imaginary part")
    refindex = epsmu.n.real
    n = qpms_emptylattice2_modes_maxfreq(&omegas_c, b1, b2, BASIS_RTOL,
	    k_, GSL_CONST_MKSA_SPEED_OF_LIGHT / refindex, maxomega)
    cdef np.ndarray[double, ndim=1] omegas = np.empty((n,), dtype=np.double)
    cdef double[::1] omegas_v = omegas
    cdef size_t i
    for i in range(n):
        omegas_v[i] = omegas_c[i]
    free(omegas_c)
    return omegas


cdef class _ScatteringSystemAtOmegaK:
    '''
    Wrapper over the C qpms_scatsys_at_omega_k_t structure
    '''
    cdef qpms_scatsys_at_omega_k_t sswk
    cdef _ScatteringSystemAtOmega ssw_pyref

    cdef qpms_scatsys_at_omega_k_t *rawpointer(self):
        return &self.sswk

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


    def ensure_finite(self):
        if self.ssw[0].ss[0].lattice_dimension != 0:
            raise NotImplementedError("Operation not supported for periodic systems")

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

    def scatter_solver(self, iri=None, k=None):
        self.check()
        cdef _ScatteringSystemAtOmegaK sswk # used only for periodic systems
        if(self.ssw[0].ss[0].lattice_dimension == 0):
            return ScatteringMatrix(self, iri=iri)
        else:
            if iri is not None:
                raise NotImplementedError("Irrep decomposition not (yet) supported for periodic systems")
            sswk = self._sswk(k)
            return ScatteringMatrix(ssw=self, sswk=sswk, iri=None)

    def _sswk(self, k):
        cdef _ScatteringSystemAtOmegaK sswk = _ScatteringSystemAtOmegaK()
        sswk.sswk.ssw = self.ssw
        sswk.sswk.k[0] = k[0]
        sswk.sswk.k[1] = k[1]
        sswk.sswk.k[2] = k[2]
        return sswk

    property fecv_size: 
        def __get__(self): return self.ss_pyref.fecv_size
    property saecv_sizes: 
        def __get__(self): return self.ss_pyref.saecv_sizes
    property irrep_names: 
        def __get__(self): return self.ss_pyref.irrep_names
    property nirreps: 
        def __get__(self): return self.ss_pyref.nirreps
    property wavenumber:
        def __get__(self): return self.ssw[0].wavenumber
        

    def modeproblem_matrix_full(self, k=None):
        self.check()
        cdef size_t flen = self.ss_pyref.s[0].fecv_size
        cdef np.ndarray[np.complex_t, ndim=2] target = np.empty(
                (flen,flen),dtype=complex, order='C')
        cdef cdouble[:,::1] target_view = target
        cdef _ScatteringSystemAtOmegaK sswk
        if k is not None:
            sswk = self._sswk(k)
            qpms_scatsyswk_build_modeproblem_matrix_full(&target_view[0][0], &sswk.sswk)
        else: 
            qpms_scatsysw_build_modeproblem_matrix_full(&target_view[0][0], self.ssw)
        return target

    def modeproblem_matrix_packed(self, qpms_iri_t iri, version='pR'):
        self.check()
        self.ensure_finite()
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

    def translation_matrix_full(self, blochvector = None):
        """Constructs full translation matrix of a scattering system at given frequency.

        Parameters
        ----------
        blochvector : array_like of shape (3,)
            Bloch vector (only for periodic systems)

        See Also
        --------
        ScatteringSystem.translation_matrix_full: translation matrix for any wavenumber

        """
        return self.ss_pyref.translation_matrix_full(wavenumber=self.wavenumber, blochvector=blochvector)

    def translation_matrix_packed(self, iri, J = QPMS_HANKEL_PLUS):
        return self.ss_pyref.translation_matrix_packed(wavenumber=self.wavenumber, iri=iri, J=J)

    def scattered_E(self, scatcoeffvector_full, evalpos, bint alt=False, btyp=QPMS_HANKEL_PLUS):
        cdef qpms_bessel_t btyp_c = BesselType(btyp)
        evalpos = np.array(evalpos, dtype=float, copy=False)
        if evalpos.shape[-1] != 3:
            raise ValueError("Last dimension of evalpos has to be 3")
        cdef np.ndarray[double,ndim=2] evalpos_a = evalpos.reshape(-1,3)
        cdef np.ndarray[dtype=complex, ndim=1] scv = np.array(scatcoeffvector_full, copy=False)
        cdef cdouble[::1] scv_view = scv
        cdef np.ndarray[complex, ndim=2] results = np.empty((evalpos_a.shape[0],3), dtype=complex)
        cdef ccart3_t res
        cdef cart3_t pos
        cdef size_t i
        for i in range(evalpos_a.shape[0]):
            pos.x = evalpos_a[i,0]
            pos.y = evalpos_a[i,1]
            pos.z = evalpos_a[i,2]
            if alt:
                res = qpms_scatsysw_scattered_E__alt(self.ssw, btyp_c, &scv_view[0], pos)
            else:
                res = qpms_scatsysw_scattered_E(self.ssw, btyp_c, &scv_view[0], pos)
            results[i,0] = res.x
            results[i,1] = res.y
            results[i,2] = res.z
        return results.reshape(evalpos.shape)


cdef class ScatteringMatrix:
    '''
    Wrapper over the C qpms_ss_LU structure that keeps the factorised mode problem matrix.
    '''
    cdef _ScatteringSystemAtOmega ssw # Here we keep the reference to the parent scattering system
    cdef _ScatteringSystemAtOmegaK sswk
    cdef qpms_ss_LU lu

    def __cinit__(self, _ScatteringSystemAtOmega ssw, sswk=None, iri=None):
        ssw.check()
        self.ssw = ssw
        if sswk is None:
            ssw.ensure_finite()
            # TODO? pre-allocate the matrix with numpy to make it transparent?
            if iri is None:
                self.lu = qpms_scatsysw_build_modeproblem_matrix_full_LU(
                        NULL, NULL, ssw.rawpointer())
            else:
                self.lu = qpms_scatsysw_build_modeproblem_matrix_irrep_packed_LU(
                        NULL, NULL, ssw.rawpointer(), iri)
        else:
            # TODO check sswk validity
            self.sswk = sswk
            self.lu = qpms_scatsyswk_build_modeproblem_matrix_full_LU(NULL, NULL, self.sswk.rawpointer())

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
        



