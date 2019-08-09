import numpy as np
from qpms_cdefs cimport *
cimport cython
import enum

# Here will be enum and dtype definitions; maybe move these to a separate file
class VSWFType(enum.IntEnum):
    ELECTRIC = QPMS_VSWF_ELECTRIC
    MAGNETIC = QPMS_VSWF_MAGNETIC
    LONGITUDINAL = QPMS_VSWF_LONGITUDINAL
    M = QPMS_VSWF_MAGNETIC
    N = QPMS_VSWF_ELECTRIC
    L = QPMS_VSWF_LONGITUDINAL

class BesselType(enum.IntEnum):
    UNDEF = QPMS_BESSEL_UNDEF
    REGULAR = QPMS_BESSEL_REGULAR
    SINGULAR = QPMS_BESSEL_SINGULAR
    HANKEL_PLUS = QPMS_HANKEL_PLUS
    HANKEL_MINUS = QPMS_HANKEL_MINUS

class PointGroupClass(enum.IntEnum):
    CN    = QPMS_PGS_CN
    S2N   = QPMS_PGS_S2N
    CNH   = QPMS_PGS_CNH
    CNV   = QPMS_PGS_CNV
    DN    = QPMS_PGS_DN
    DND   = QPMS_PGS_DND
    DNH   = QPMS_PGS_DNH
    T     = QPMS_PGS_T
    TD    = QPMS_PGS_TD
    TH    = QPMS_PGS_TH
    O     = QPMS_PGS_O
    OH    = QPMS_PGS_OH
    I     = QPMS_PGS_I
    IH    = QPMS_PGS_IH
    CINF  = QPMS_PGS_CINF
    CINFH = QPMS_PGS_CINFH
    CINFV = QPMS_PGS_CINFV
    DINF  = QPMS_PGS_DINF
    DINFH = QPMS_PGS_DINFH
    SO3   = QPMS_PGS_SO3
    O3    = QPMS_PGS_O3

try:
    class DebugFlags(enum.IntFlag): # Should be IntFlag if python version >= 3.6
        MISC = QPMS_DBGMSG_MISC
        THREADS = QPMS_DBGMSG_THREADS
    has_IntFlag = True
except AttributeError: # For old versions of enum, use IntEnum instead
    class DebugFlags(enum.IntEnum): 
        MISC = QPMS_DBGMSG_MISC
        THREADS = QPMS_DBGMSG_THREADS
    has_IntFlag = False

def dbgmsg_enable(qpms_dbgmsg_flags types):
    flags = qpms_dbgmsg_enable(types)
    return DebugFlags(flags) if has_IntFlag else flags
def dbgmsg_disable(qpms_dbgmsg_flags types):
    flags = qpms_dbgmsg_disable(types)
    return DebugFlags(flags) if has_IntFlag else flags
def dbgmsg_active():
    flags = qpms_dbgmsg_enable(<qpms_dbgmsg_flags>0)
    return DebugFlags(flags) if has_IntFlag else flags

#import re # TODO for crep methods?

#cimport openmp
#openmp.omp_set_dynamic(1)

## Auxillary function for retrieving the "meshgrid-like" indices; inc. nmax
@cython.boundscheck(False)
def get_mn_y(int nmax):
    """
    Auxillary function for retreiving the 'meshgrid-like' indices from the flat indexing; 
    inc. nmax.
    ('y to mn' conversion)
    
    Parameters
    ----------

    nmax : int
        The maximum order to which the VSWFs / Legendre functions etc. will be evaluated.
        
    Returns
    -------
    
    output : (m, n)
        Tuple of two arrays of type np.array(shape=(nmax*nmax + 2*nmax), dtype=np.int),
        where [(m[y],n[y]) for y in range(nmax*nmax + 2*nma)] covers all possible 
        integer pairs n >= 1, -n <= m <= n.
    """
    cdef Py_ssize_t nelems = nmax * nmax + 2 * nmax
    cdef np.ndarray[np.int_t,ndim=1] m_arr = np.empty([nelems], dtype=np.int)
    cdef np.ndarray[np.int_t,ndim=1] n_arr = np.empty([nelems], dtype=np.int)
    cdef Py_ssize_t i = 0
    cdef np.int_t n, m
    for n in range(1,nmax+1):
        for m in range(-n,n+1):
            m_arr[i] = m
            n_arr[i] = n
            i = i + 1
    return (m_arr, n_arr)

def get_nelem(unsigned int lMax):
    return lMax * (lMax + 2)

def get_y_mn_unsigned(int nmax): 
    """
    Auxillary function for mapping 'unsigned m', n indices to the flat y-indexing.
    For use with functions as scipy.special.lpmn, which have to be evaluated separately
    for positive and negative m.
    
    Parameters
    ----------

    nmax : int
        The maximum order to which the VSWFs / Legendre functions etc. will be evaluated.
        
    output : (ymn_plus, ymn_minus)
        Tuple of two arrays of shape (nmax+1,nmax+1), containing the flat y-indices corresponding
        to the respective (m,n) and (-m,n). The elements for which |m| > n are set to -1.
        (Therefore, the caller must not use those elements equal to -1.)
    """
    cdef np.ndarray[np.intp_t, ndim=2] ymn_plus = np.full((nmax+1,nmax+1),-1, dtype=np.intp)
    cdef np.ndarray[np.intp_t, ndim=2] ymn_minus = np.full((nmax+1,nmax+1),-1, dtype=np.intp)
    cdef Py_ssize_t i = 0
    cdef np.int_t n, m
    for n in range(1,nmax+1):
        for m in range(-n,0):
            ymn_minus[-m,n] = i
            i = i + 1
        for m in range(0,n+1):
            ymn_plus[m,n] = i
            i = i + 1
    return(ymn_plus, ymn_minus)


def tlm2uvswfi(t, l, m):
    ''' TODO doc
    And TODO this should rather be an ufunc.
    '''
    # Very low-priority TODO: add some types / cythonize
    if isinstance(t, int) and isinstance(l, int) and isinstance(m, int):
        return qpms_tmn2uvswfi(t, m, l)
    elif len(t) == len(l) and len(t) == len(m):
        u = list()
        for i in range(len(t)):
            if not (t[i] % 1 == 0 and l[i] % 1 == 0 and m[i] % 1 == 0): # maybe not the best check possible, though
                raise ValueError # TODO error message
            u.append(qpms_tmn2uvswfi(t[i],m[i],l[i]))
        return u
    else:
        print(len(t), len(l), len(m))
        raise ValueError("Lengths of the t,l,m arrays must be equal, but they are %d, %d, %d."
                % (len(t), len(l), len(m)))


def uvswfi2tlm(u):
    ''' TODO doc
    and TODO this should rather be an ufunc.
    '''
    cdef qpms_vswf_type_t t
    cdef qpms_l_t l
    cdef qpms_m_t m
    cdef size_t i
    if isinstance(u, (int, np.ulonglong)):
        if (qpms_uvswfi2tmn(u, &t, &m, &l) != QPMS_SUCCESS):
            raise ValueError("Invalid uvswf index")
        return (t, l, m)
    else:
        ta = list()
        la = list()
        ma = list()
        for i in range(len(u)):
            if (qpms_uvswfi2tmn(u[i], &t, &m, &l) != QPMS_SUCCESS):
                raise ValueError("Invalid uvswf index")
            ta.append(t)
            la.append(l)
            ma.append(m)
        return (ta, la, ma)

