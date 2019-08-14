import numpy as np
import enum
from .cycommon import get_mn_y, tlm2uvswfi

class VSWFNorm(enum.IntEnum):
    # TODO try to make this an enum.IntFlag if supported
    # TODO add the other flags from qpms_normalisation_t as well
    UNNORM = QPMS_NORMALISATION_NORM_NONE
    UNNORM_CS = QPMS_NORMALISATION_NORM_NONE | QPMS_NORMALISATION_CSPHASE
    POWERNORM = QPMS_NORMALISATION_NORM_POWER
    POWERNORM_CS = QPMS_NORMALISATION_NORM_POWER | QPMS_NORMALISATION_CSPHASE
    SPHARMNORM = QPMS_NORMALISATION_NORM_SPHARM
    SPHARMNORM_CS = QPMS_NORMALISATION_NORM_SPHARM | QPMS_NORMALISATION_CSPHASE
    UNDEF = QPMS_NORMALISATION_UNDEF

cdef class BaseSpec:
    '''Cython wrapper over qpms_vswf_set_spec_t.

    It should be kept immutable. The memory is managed by numpy/cython, not directly by the C functions, therefore
    whenever used in other wrapper classes that need the pointer
    to qpms_vswf_set_spec_t, remember to set a (private, probably immutable) reference to qpms.basespec to ensure
    correct reference counting and garbage collection.
    '''
    #cdef qpms_vswf_set_spec_t s # in pxd
    #cdef np.ndarray __ilist     # in pxd
    #cdef const qpms_uvswfi_t[:] __ilist

    def __cinit__(self, *args, **kwargs):
        cdef const qpms_uvswfi_t[:] ilist_memview
        if len(args) == 0:
            if 'lMax' in kwargs.keys(): # if only lMax is specified, create the 'usual' definition in ('E','M') order
                lMax = kwargs['lMax']
                my, ny = get_mn_y(lMax)
                nelem = len(my)
                tlist = nelem * (QPMS_VSWF_ELECTRIC,) + nelem * (QPMS_VSWF_MAGNETIC,)
                mlist = 2*list(my)
                llist = 2*list(ny)
                ilist = tlm2uvswfi(tlist,llist,mlist)
            else:
                raise ValueError
        else: # len(args) > 0:
            ilist = args[0]
            #self.__ilist = np.array(args[0], dtype=qpms_uvswfi_t, order='C', copy=True) # FIXME define the dtypes at qpms_cdef.pxd level
        self.__ilist = np.array(ilist, dtype=np.ulonglong, order='C', copy=True)
        self.__ilist.setflags(write=False)
        ilist_memview = self.__ilist
        self.s.ilist = &ilist_memview[0]
        self.s.n = len(self.__ilist)
        self.s.capacity = 0 # is this the best way?
        if 'norm' in kwargs.keys():
            self.s.norm = kwargs['norm']
        else:
            self.s.norm = <qpms_normalisation_t>(QPMS_NORMALISATION_NORM_POWER | QPMS_NORMALISATION_CSPHASE)
        # set the other metadata
        cdef qpms_l_t l
        self.s.lMax_L = -1
        cdef qpms_m_t m
        cdef qpms_vswf_type_t t
        for i in range(self.s.n):
            if(qpms_uvswfi2tmn(ilist_memview[i], &t, &m, &l) != QPMS_SUCCESS):
                raise ValueError("Invalid uvswf index")
            if (t == QPMS_VSWF_ELECTRIC):
                self.s.lMax_N = max(self.s.lMax_N, l)
            elif (t == QPMS_VSWF_MAGNETIC):
                self.s.lMax_M = max(self.s.lMax_M, l)
            elif (t == QPMS_VSWF_LONGITUDINAL):
                self.s.lMax_L = max(self.s.lMax_L, l)
            else:
                raise ValueError # If this happens, it's probably a bug, as it should have failed already at qpms_uvswfi2tmn
            self.s.lMax = max(self.s.lMax, l)

    def tlm(self):
        cdef const qpms_uvswfi_t[:] ilist_memview = <qpms_uvswfi_t[:self.s.n]> self.s.ilist
        #cdef qpms_vswf_type_t[:] t = np.empty(shape=(self.s.n,), dtype=qpms_vswf_type_t) # does not work, workaround:
        cdef size_t i
        cdef np.ndarray ta = np.empty(shape=(self.s.n,), dtype=np.intc)
        cdef int[:] t = ta 
        #cdef qpms_l_t[:] l = np.empty(shape=(self.s.n,), dtype=qpms_l_t) # FIXME explicit dtype again
        cdef np.ndarray la = np.empty(shape=(self.s.n,), dtype=np.intc) 
        cdef qpms_l_t[:] l = la 
        #cdef qpms_m_t[:] m = np.empty(shape=(self.s.n,), dtype=qpms_m_t) # FIXME explicit dtype again
        cdef np.ndarray ma =  np.empty(shape=(self.s.n,), dtype=np.intc) 
        cdef qpms_m_t[:] m = ma
        for i in range(self.s.n):
            qpms_uvswfi2tmn(self.s.ilist[i], <qpms_vswf_type_t*>&t[i], &m[i], &l[i])
        return (ta, la, ma)

    def m(self): # ugly
        return self.tlm()[2]

    def t(self): # ugly
        return self.tlm()[0]

    def l(self): # ugly
        return self.tlm()[1]

    def __len__(self):
        return self.s.n

    def __getitem__(self, key):
        # TODO raise correct errors (TypeError on bad type of key, IndexError on exceeding index)
        return self.__ilist[key]

    property ilist:
        def __get__(self):
            return self.__ilist

    cdef qpms_vswf_set_spec_t *rawpointer(BaseSpec self):
        '''Pointer to the qpms_vswf_set_spec_t structure.
        Don't forget to reference the BaseSpec object itself when storing the pointer anywhere!!!
        '''
        return &(self.s)

    property rawpointer:
        def __get__(self):
            return <uintptr_t> &(self.s)
    
    property norm:
        def __get__(self):
            return VSWFNorm(self.s.norm)

default_bspec = BaseSpec(lMax=2)
