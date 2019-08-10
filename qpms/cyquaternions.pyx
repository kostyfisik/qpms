from .cybspec cimport BaseSpec
from .qpms_cdefs cimport *
import cmath
import math


def complex_crep(complex c, parentheses = False, shortI = True, has_Imaginary = False):
    '''
    Return a C-code compatible string representation of a (python) complex number.
    '''
    return ( ('(' if parentheses else '')
            + repr(c.real)
            + ('+' if math.copysign(1, c.imag) >= 0 else '')
            + repr(c.imag)
            + ('*I' if shortI else '*_Imaginary_I' if has_Imaginary else '*_Complex_I')
            + (')' if parentheses else '')
        )


cdef class CQuat:
    '''
    Wrapper of the qpms_quat_t object, with the functionality
    to evaluate Wigner D-matrix elements.
    '''
    # cdef readonly qpms_quat_t q # pxd

    def __cinit__(self, double w, double x, double y, double z):
        cdef qpms_quat4d_t p
        p.c1 = w
        p.ci = x
        p.cj = y
        p.ck = z
        self.q = qpms_quat_2c_from_4d(p)

    def copy(self):
        res = CQuat(0,0,0,0)
        res.q = self.q
        return res

    def __repr__(self): # TODO make this look like a quaternion with i,j,k
        return repr(self.r)

    def __add__(CQuat self, CQuat other):
        # TODO add real numbers
        res = CQuat(0,0,0,0)
        res.q = qpms_quat_add(self.q, other.q)
        return res

    def __mul__(self, other):
        res = CQuat(0,0,0,0)
        if isinstance(self, CQuat):
            if isinstance(other, CQuat):
                res.q = qpms_quat_mult(self.q, other.q)
            elif isinstance(other, (int, float)):
                res.q = qpms_quat_rscale(other, self.q)
            else: return NotImplemented
        elif isinstance(self, (int, float)):
            if isinstance(other, CQuat):
                res.q = qpms_quat_rscale(self, other.q)
            else: return NotImplemented
        return res

    def __neg__(CQuat self):
        res = CQuat(0,0,0,0)
        res.q = qpms_quat_rscale(-1, self.q)
        return res

    def __sub__(CQuat self, CQuat other):
        res = CQuat(0,0,0,0)
        res.q = qpms_quat_add(self.q, qpms_quat_rscale(-1,other.q))
        return res

    def __abs__(self):
        return qpms_quat_norm(self.q)

    def norm(self):
        return qpms_quat_norm(self.q)

    def imnorm(self):
        return qpms_quat_imnorm(self.q)

    def exp(self):
        res = CQuat(0,0,0,0)
        res.q = qpms_quat_exp(self.q)
        return res

    def log(self):
        res = CQuat(0,0,0,0)
        res.q = qpms_quat_exp(self.q)
        return res

    def __pow__(CQuat self, double other, _):
        res = CQuat(0,0,0,0)
        res.q = qpms_quat_pow(self.q, other)
        return res

    def normalise(self):
        res = CQuat(0,0,0,0)
        res.q = qpms_quat_normalise(self.q)
        return res

    def isclose(CQuat self, CQuat other, rtol=1e-5, atol=1e-8):
        '''
        Checks whether two quaternions are "almost equal".
        '''
        return abs(self - other) <= (atol + rtol * abs(other))

    property c:
        '''
        Quaternion representation as two complex numbers
        '''
        def __get__(self):
            return (self.q.a, self.q.b)
        def __set__(self, RaRb):
            self.q.a = RaRb[0]
            self.q.b = RaRb[1]

    property r:
        '''
        Quaternion representation as four real numbers
        '''
        def __get__(self):
            cdef qpms_quat4d_t p
            p = qpms_quat_4d_from_2c(self.q)
            return (p.c1, p.ci, p.cj, p.ck)
        def __set__(self, wxyz):
            cdef qpms_quat4d_t p
            p.c1 = wxyz[0]
            p.ci = wxyz[1]
            p.cj = wxyz[2]
            p.ck = wxyz[3]
            self.q = qpms_quat_2c_from_4d(p)

    def crepr(self):
        '''
        Returns a string that can be used in C code to initialise a qpms_irot3_t
        '''
        return '{' + complex_crep(self.q.a) + ', ' + complex_crep(self.q.b)  + '}'

    def wignerDelem(self, qpms_l_t l, qpms_m_t mp, qpms_m_t m):
        '''
        Returns an element of a bosonic Wigner matrix.
        '''
        # don't crash on bad l, m here
        if (abs(m) > l or abs(mp) > l):
            return 0
        return qpms_wignerD_elem(self.q, l, mp, m)
    
    @staticmethod
    def from_rotvector(vec):
        if vec.shape != (3,):
            raise ValueError("Single 3d vector expected")
        res = CQuat()
        cdef cart3_t v
        v.x = vec[0]
        v.y = vec[1]
        v.z = vec[2]
        res.q = qpms_quat_from_rotvector(v)
        return res

cdef class IRot3:
    '''
    Wrapper over the C type qpms_irot3_t.
    '''
    #cdef readonly qpms_irot3_t qd

    def __cinit__(self, *args): 
        '''
        TODO doc
        '''
        # TODO implement a constructor with
        #  - tuple as argument ...?
        if (len(args) == 0): # no args, return identity
            self.qd.rot.a = 1
            self.qd.rot.b = 0
            self.qd.det = 1
        elif (len(args) == 2 and isinstance(args[0], CQuat) and isinstance(args[1], (int, float))):
            # The original __cinit__(self, CQuat q, short det) constructor
            q = args[0]
            det = args[1]
            if (det != 1 and det != -1):
                raise ValueError("Improper rotation determinant has to be 1 or -1")
            self.qd.rot = q.normalise().q
            self.qd.det = det
        elif (len(args) == 1 and isinstance(args[0], IRot3)):
            # Copy
            self.qd = args[0].qd
        elif (len(args) == 1 and isinstance(args[0], CQuat)):
            # proper rotation from a quaternion
            q = args[0]
            det = 1
            self.qd.rot = q.normalise().q
            self.qd.det = det
        else:
            raise ValueError('Unsupported constructor arguments')

    cdef void cset(self, qpms_irot3_t qd):
        self.qd = qd

    def copy(self):
        res = IRot3(CQuat(1,0,0,0),1)
        res.qd = self.qd
        return res

    property rot:
        '''
        The proper rotation part of the IRot3 type.
        '''
        def __get__(self):
            res = CQuat(0,0,0,0)
            res.q = self.qd.rot
            return res
        def __set__(self, CQuat r):
            # TODO check for non-zeroness and throw an exception if norm is zero
            self.qd.rot = r.normalise().q

    property det:
        '''
        The determinant of the improper rotation.
        '''
        def __get__(self):
            return self.qd.det
        def __set__(self, d):
            d = int(d)
            if (d != 1 and d != -1):
                raise ValueError("Improper rotation determinant has to be 1 or -1")
            self.qd.det = d

    def __repr__(self): # TODO make this look like a quaternion with i,j,k
        return '(' + repr(self.rot) + ', ' + repr(self.det) + ')'

    def crepr(self):
        '''
        Returns a string that can be used in C code to initialise a qpms_irot3_t
        '''
        return '{' + self.rot.crepr() + ', ' + repr(self.det) + '}'

    def __mul__(IRot3 self, IRot3 other):
        res = IRot3(CQuat(1,0,0,0), 1) 
        res.qd = qpms_irot3_mult(self.qd, other.qd)
        return res

    def __pow__(IRot3 self, n, _):
        cdef int nint
        if (n % 1 == 0):
            nint = n
        else:
            raise ValueError("The exponent of an IRot3 has to have an integer value.")
        res = IRot3(CQuat(1,0,0,0), 1)
        res.qd = qpms_irot3_pow(self.qd, n)
        return res

    def isclose(IRot3 self, IRot3 other, rtol=1e-5, atol=1e-8):
        '''
        Checks whether two (improper) rotations are "almost equal".
        Returns always False if the determinants are different.
        '''
        if self.det != other.det: 
            return False
        return (self.rot.isclose(other.rot, rtol=rtol, atol=atol)
                # unit quaternions are a double cover of SO(3), i.e.
                # minus the same quaternion represents the same rotation
                or self.rot.isclose(-(other.rot), rtol=rtol, atol=atol)
            )

    # Several 'named constructors' for convenience
    @staticmethod
    def inversion():
        '''
        Returns an IRot3 object representing the 3D spatial inversion.
        '''
        r = IRot3()
        r.det = -1
        return r

    @staticmethod
    def zflip():
        '''
        Returns an IRot3 object representing the 3D xy-plane mirror symmetry (z axis sign flip).
        '''
        r = IRot3()
        r.rot = CQuat(0,0,0,1) # π-rotation around z-axis
        r.det = -1 # inversion
        return r

    @staticmethod
    def yflip():
        '''
        Returns an IRot3 object representing the 3D xz-plane mirror symmetry (y axis sign flip).
        '''
        r = IRot3()
        r.rot = CQuat(0,0,1,0) # π-rotation around y-axis
        r.det = -1 # inversion
        return r

    @staticmethod
    def xflip():
        '''
        Returns an IRot3 object representing the 3D yz-plane mirror symmetry (x axis sign flip).
        '''
        r = IRot3()
        r.rot = CQuat(0,1,0,0) # π-rotation around x-axis
        r.det = -1 # inversion
        return r

    @staticmethod
    def zrotN(int n):
        '''
        Returns an IRot3 object representing a \f$ C_n $\f rotation (around the z-axis).
        '''
        r = IRot3()
        r.rot = CQuat(math.cos(math.pi/n),0,0,math.sin(math.pi/n))
        return r

    @staticmethod
    def identity():
        '''
        An alias for the constructor without arguments; returns identity.
        '''
        return IRot3()

    def as_uvswf_matrix(IRot3 self, BaseSpec bspec):
        '''
        Returns the uvswf representation of the current transform as a numpy array
        '''
        cdef ssize_t sz = len(bspec)
        cdef np.ndarray m = np.empty((sz, sz), dtype=complex, order='C') # FIXME explicit dtype
        cdef cdouble[:, ::1] view = m
        qpms_irot3_uvswfi_dense(&view[0,0], bspec.rawpointer(), self.qd)
        return m

