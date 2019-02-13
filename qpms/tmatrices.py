import numpy as np
import quaternion, spherical_functions as sf # because of the Wigner matrices. These imports are SLOW.
import re
from scipy import interpolate
from scipy.constants import hbar, e as eV, pi, c
from qpms_c import get_mn_y, get_nelem
ň = np.newaxis
from .types import NormalizationT, TMatrixSpec

# Transformations of spherical bases
def WignerD_mm(l, quat):
    """
    Calculates Wigner D matrix (as an numpy (2*l+1,2*l+1)-shaped array) 
    for order l, and a rotation given by quaternion quat.
    
    This represents the rotation of spherical vector basis
    TODO doc
    """
    
    indices = np.array([ [l,i,j] for i in range(-l,l+1) for j in range(-l,l+1)])
    Delems = sf.Wigner_D_element(quat, indices).reshape(2*l+1,2*l+1)
    return Delems

def WignerD_mm_fromvector(l, vect):
    """
    TODO doc
    """
    return WignerD_mm(l, quaternion.from_rotation_vector(vect))


def WignerD_yy(lmax, quat):
    """
    TODO doc
    """
    my, ny = get_mn_y(lmax)
    Delems = np.zeros((len(my),len(my)),dtype=complex)
    b_in = 0
    e_in =  None
    for l in range(1,lmax+1):
        e_in = b_in + 2*l+1
        Delems[b_in:e_in,b_in:e_in] = WignerD_mm(l, quat)
        b_in = e_in
    return Delems
    
def WignerD_yy_fromvector(lmax, vect):
    """
    TODO doc
    """
    return WignerD_yy(lmax, quaternion.from_rotation_vector(vect))

def xflip_yy(lmax):
    """
    TODO doc
    xflip = δ(m + m') δ(l - l')  
    (i.e. ones on the (m' m) antidiagonal
    """
    my, ny = get_mn_y(lmax)
    elems = np.zeros((len(my),len(my)),dtype=int)
    b_in = 0
    e_in =  None
    for l in range(1,lmax+1):
        e_in = b_in + 2*l+1
        elems[b_in:e_in,b_in:e_in] = np.eye(2*l+1)[::-1,:]
        b_in = e_in
    return elems

def xflip_tyy(lmax):
    fl_yy = xflip_yy(lmax)
    return np.array([fl_yy,-fl_yy])

def xflip_tyty(lmax):
    fl_yy = xflip_yy(lmax)
    nelem = fl_yy.shape[0]
    fl_tyty = np.zeros((2,nelem,2,nelem),dtype=int)
    fl_tyty[0,:,0,:] = fl_yy
    fl_tyty[1,:,1,:] = -fl_yy
    return fl_tyty

def yflip_yy(lmax):
    """
    TODO doc
    yflip = rot(z,pi/2) * xflip * rot(z,-pi/2)
          = δ(m + m') δ(l - l') * (-1)**m
    """
    my, ny = get_mn_y(lmax)
    elems = xflip_yy(lmax)
    elems[(my % 2)==1] = elems[(my % 2)==1] * -1 # Obvious sign of tiredness (this is correct but ugly; FIXME)
    return elems

def yflip_tyy(lmax):
    fl_yy = yflip_yy(lmax)
    return np.array([fl_yy,-fl_yy])

def yflip_tyty(lmax):
    fl_yy = yflip_yy(lmax)
    nelem = fl_yy.shape[0]
    fl_tyty = np.zeros((2,nelem,2,nelem),dtype=int)
    fl_tyty[0,:,0,:] = fl_yy
    fl_tyty[1,:,1,:] = -fl_yy
    return fl_tyty

def zflip_yy(lmax):
    """
    TODO doc
    zflip = (-1)^(l+m)
    """
    my, ny = get_mn_y(lmax)
    elems = np.zeros((len(my), len(my)), dtype=int)
    b_in = 0
    e_in = None
    for l in range(1,lmax+1):
        e_in = b_in + 2*l+1
        elems[b_in:e_in,b_in:e_in] = np.diag([(-1)**i for i in range(e_in-b_in)])
        b_in = e_in
    return elems

def zflip_tyy(lmax):
    fl_yy = zflip_yy(lmax)
    return np.array([fl_yy,-fl_yy])

def zflip_tyty(lmax):
    fl_yy = zflip_yy(lmax)
    nelem = fl_yy.shape[0]
    fl_tyty = np.zeros((2,nelem,2,nelem),dtype=int)
    fl_tyty[0,:,0,:] = fl_yy
    fl_tyty[1,:,1,:] = -fl_yy
    return fl_tyty

def zrotN_yy(N, lMax):
    return WignerD_yy_fromvector(lMax, np.array([0,0,pi * (2/N)]))

def op_yy2tyty(yyop):
    '''
    Broadcasts an yy operator to tyty operator without considering mirroring effects.
    Good (maybe only) for rotations.
    '''
    return np.moveaxis(np.eye(2)[:,:,ň,ň] * yyop, 2,1)

def zrotN_tyty(N, lMax):
    return op_yy2tyty(zrotN_yy(N, lMax))

def parity_yy(lmax):
    """
    Parity operator (flip in x,y,z)
    parity = (-1)**l
    """
    my, ny = get_mn_y(lmax)
    return np.diag((-1)**ny)

# BTW parity (xyz-flip) is simply (-1)**ny



#----------------------------------------------------#
# Loading T-matrices from scuff-tmatrix output files #
#----------------------------------------------------#

# We don't really need this particular function anymore, but...
def _scuffTMatrixConvert_EM_01(EM):
    #print(EM)
    if (EM == b'E'):
        return 1
    elif (EM == b'M'):
        return 0
    else:
        return None

def loadScuffTMatrices(fileName, normalisation = 1, version = 'old', freqscale = None, order = None):
    """
    TODO doc

    version describes version of scuff-em. It is either 'old' or 'new'. 
    default order is ('N','M') with 'old' version, ('M','N') with 'new'
    """
    oldversion = (version == 'old')
    μm = 1e-6
    table = np.genfromtxt(fileName, 
                  converters={1: _scuffTMatrixConvert_EM_01, 4: _scuffTMatrixConvert_EM_01} if oldversion else None,
                  skip_header = 0 if oldversion else 5,
                  usecols = None if oldversion else (0, 2, 3, 4, 6, 7, 8, 9, 10),
                  dtype=[('freq', '<f8'), 
                         ('outc_type', '<i8'), ('outc_l', '<i8'), ('outc_m', '<i8'), 
                         ('inc_type', '<i8'), ('inc_l', '<i8'), ('inc_m', '<i8'), 
                         ('Treal', '<f8'), ('Timag', '<f8')] 
                   if oldversion else 
                        [('freq', '<f8'), 
                         ('outc_l', '<i8'), ('outc_m', '<i8'),  ('outc_type', '<i8'),
                         ('inc_l', '<i8'), ('inc_m', '<i8'), ('inc_type', '<i8'), 
                         ('Treal', '<f8'), ('Timag', '<f8')]
                    )
    lMax=np.max(table['outc_l'])
    my,ny = get_mn_y(lMax)
    nelem = len(ny)
    TMatrix_sz = nelem**2 * 4 # number of rows for each frequency: nelem * nelem spherical incides, 2 * 2 E/M types
    freqs_weirdunits = table['freq'][::TMatrix_sz].copy()
    freqs = freqs_weirdunits * c / μm
    # The iteration in the TMatrix file goes in this order (the last one iterates fastest, i.e. in the innermost loop):
    # freq outc_l outc_m outc_type inc_l inc_m inc_type
    # The l,m mapping is the same as is given by my get_mn_y function, so no need to touch that
    TMatrices_tmp_real = table['Treal'].reshape(len(freqs), nelem, 2, nelem, 2)
    TMatrices_tmp_imag = table['Timag'].reshape(len(freqs), nelem, 2, nelem, 2)
    # There are two přoblems with the previous matrices. First, we want to have the
    # type indices first, so we want a shape (len(freqs), 2, nelem, 2, nelem) as in the older code.
    # Second, M-waves come first, so they have now 0-valued index, and E-waves have 1-valued index,
    # which we want to be inverted.
    TMatrices = np.zeros((len(freqs),2,nelem,2,nelem),dtype=complex)
    reorder = (0,1)
    if ((order == ('N', 'M')) and not oldversion): # reverse order for the new version
        reorder = (1,0) 
    # TODO reverse order for the old version...
    for inc_type in reorder:
        for outc_type in reorder:
            TMatrices[:,outc_type,:,inc_type,:] = TMatrices_tmp_real[:,:,outc_type,:,inc_type]+1j*TMatrices_tmp_imag[:,:,outc_type,:,inc_type]
    # IMPORTANT: now we are going from Reid's/Kristensson's/Jackson's/whoseever convention to Taylor's convention
    # TODO make these consistent with what is defined in qpms_types.h (implement all possibilities)
    if normalisation == 1:
        TMatrices[:,:,:,:,:] = TMatrices[:,:,:,:,:] * np.sqrt(ny*(ny+1))[ň,ň,ň,ň,:] / np.sqrt(ny*(ny+1))[ň,ň,:,ň,ň]
    elif normalisation == 2: # Kristensson?
        pass
    if freqscale is not None:
        freqs *= freqscale
        freqs_weirdunits *= freqscale
    return (TMatrices, freqs, freqs_weirdunits, lMax)


# misc tensor maniputalion
def apply_matrix_left(matrix, tensor, axis):
    """
    TODO doc
    Apply square matrix to a given axis of a tensor, so that the result retains the shape
    of the original tensor. The summation goes over the second index of the matrix and the
    given tensor axis.
    """
    tmp = np.tensordot(matrix, tensor, axes=(-1,axis))
    return np.moveaxis(tmp, 0, axis)

def apply_ndmatrix_left(matrix,tensor,axes):
    """
    Generalized apply_matrix_left, the matrix can have more (2N) abstract dimensions,
    like M[i,j,k,...z,i,j,k,...,z]. N axes have to be specified in a tuple, corresponding
    to the axes 0,1,...N-1 of the matrix
    """
    N = len(axes)
    matrix = np.tensordot(matrix, tensor, axes=([-N+axn for axn in range(N)],axes))
    matrix = np.moveaxis(matrix, range(N), axes)
    return matrix

def apply_ndmatrix_right(tensor, matrix, axes):
    """
    Right-side analogue of apply_ndmatrix_lift.
    Multiplies a tensor with a 2N-dimensional matrix, conserving the axis order.
    """
    N = len(axes)
    matrix = np.tensordot(tensor, matrix, axes = (axes, range(N)))
    matrix = np.moveaxis(matrix, [-N+axn for axn in range(N)], axes)
    return matrix

def ndmatrix_Hconj(matrix):
    """
    For 2N-dimensional matrix, swap the first N and last N matrix, and complex conjugate.
    """
    twoN = len(matrix.shape)
    if not twoN % 2 == 0:
        raise ValueError("The matrix has to have even number of axes.")
    N = twoN//2
    matrix = np.moveaxis(matrix, range(N), range(N, 2*N))
    return matrix.conj()

def symz_indexarrays(lMax, npart = 1):
    """
    Returns indices that are used for separating the in-plane E ('TE' in the photonic crystal
    jargon) and perpendicular E ('TM' in the photonic crystal jargon) modes
    in the z-mirror symmetric systems.

    Parameters
    ----------
    lMax : int
        The maximum degree cutoff for the T-matrix to which these indices will be applied.

    npart : int
        Number of particles (TODO better description)

    Returns
    -------
    TEč, TMč : (npart * 2 * nelem)-shaped bool ndarray
        Mask arrays corresponding to the 'TE' and 'TM' modes, respectively.
    """
    my, ny = get_mn_y(lMax)
    nelem = len(my)
    ž = np.arange(2*nelem) # single particle spherical wave indices
    tž = ž // nelem # tž == 0: electric waves, tž == 1: magnetic waves
    mž = my[ž%nelem]
    nž = ny[ž%nelem]
    TEž = ž[(mž+nž+tž) % 2 == 0]
    TMž = ž[(mž+nž+tž) % 2 == 1]

    č = np.arange(npart*2*nelem) # spherical wave indices for multiple particles (e.g. in a unit cell)
    žč = č % (2* nelem)
    tč = tž[žč]
    mč = mž[žč]
    nč = nž[žč]
    TEč = č[(mč+nč+tč) % 2 == 0]
    TMč = č[(mč+nč+tč) % 2 == 1]
    return (TEč, TMč)


"""
Processing T-matrix related operations from scripts
===================================================

see also scripts_common.py

The unit cell is defined by a dict particle_specs and a list TMatrix_specs.
This particular module has to provide the T-matrices according to what is defined
in TMatrix_specs

TMatrix_specs is a list of tuples (lMax_override, TMatrix_path, ops)
where 
  - TMatrix_path is path to the file generated by scuff-tmatrix
  - lMax_override is int or None; if it is int and less than the lMax found from the T-matrix file,
    lMax_override is used as the order cutoff for the output T-matrix.
  - ops is an iterable of tuples (optype, opargs) where currently optype can be 'sym' or 'tr'
    for symmetrization operation or some other transformation.
"""

#TODO FEATURE more basic group symmetry operations, cf. http://symmetry.otterbein.edu/tutorial/index.html

# This is for finite „fractional“ rotations along the z-axis (mCN means rotation of 2π*(m/N))
reCN = re.compile('(\d*)C(\d+)') # TODO STYLE make this regexp also accept the 3*C_5-type input, eqiv. to 3C5

def get_TMatrix_fromspec(tmatrix_spec):
    ''' TODO doc 
    returns (TMatrices, freqs, lMax)
    '''
    lMax_override, tmpath, ops = tmatrix_spec
    TMatrices, freqs, freqs_weirdunits, lMax = loadScuffTMatrices(tmpath)
    if lMax_override is not None and (lMax_override < lMax_orig):
        nelem = get_nelem(lMax_override)
        TMatrices = TMatrices[...,0:nelem,:,0:nelem]
        lMax = lMax_override
    
    for (optype, opargs) in ops:
        if optype == 'sym':
            mCN = reCN.match(opargs)
            if opargs == 'C2' or opargs == 'C_2':
                opmat = apply_matrix_left(yflip_yy(lMax), xflip_yy(lMax), -1)
                TMatrices = (TMatrices + apply_matrix_left(opmat, apply_matrix_left(opmat, TMatrices, -3), -1))/2
            elif opargs == 'σ_x' or opargs == 'σx':
                opmat = xflip_tyty(lMax)
                TMatrices = (TMatrices + apply_ndmatrix_left(opmat, apply_ndmatrix_left(opmat, TMatrices, (-4,-3)),(-2,-1)))/2
            elif opargs == 'σ_y' or opargs == 'σy':
                opmat = yflip_tyty(lMax)
                TMatrices = (TMatrices + apply_ndmatrix_left(opmat, apply_ndmatrix_left(opmat, TMatrices, (-4,-3)),(-2,-1)))/2
            elif opargs == 'σ_z' or opargs == 'σz':
                opmat = zflip_tyty(lMax)
                TMatrices = (TMatrices + apply_ndmatrix_left(opmat, apply_ndmatrix_left(opmat, TMatrices, (-4,-3)),(-2,-1)))/2
            elif mCN:
                rotN = int(mCN.group(2)) # the possible m is ignored
                TMatrix_contribs = np.empty((rotN,)+TMatrices.shape, dtype=np.complex_)
                for i in range(rotN):
                    rotangle = 2*np.pi*i / rotN
                    rot = WignerD_yy_fromvector(lMax, np.array([0,0,rotangle]))
                    rotinv = WignerD_yy_fromvector(lMax, np.array([0,0,-rotangle]))
                    TMatrix_contribs[i] = apply_matrix_left(rot, apply_matrix_left(rotinv, TMatrices, -3), -1)
                TMatrices = np.sum(TMatrix_contribs, axis=0) / rotN
            elif opargs == 'C_inf' or opargs == 'Cinf' or opargs == 'C_∞' or opargs == 'C∞':
                raise ValueError('not implemented: ', opargs)  # TODO FEATURE
            else:
                raise ValueError('not implemented: ', opargs)
        elif optype == 'tr':
            mCN = reCN.match(opargs)
            if opargs == 'C2' or opargs == 'C_2':
                opmat = apply_matrix_left(yflip_yy(lMax), xflip_yy(lMax), -1)
                TMatrices = apply_matrix_left(opmat, apply_matrix_left(opmat, TMatrices, -3), -1)/2
            elif opargs == 'σ_x' or opargs == 'σx':
                opmat = xflip_tyty(lMax)
                TMatrices = apply_ndmatrix_left(opmat, apply_ndmatrix_left(opmat, TMatrices, (-4,-3)),(-2,-1))/2
            elif opargs == 'σ_y' or opargs == 'σy':
                opmat = yflip_tyty(lMax)
                TMatrices = apply_ndmatrix_left(opmat, apply_ndmatrix_left(opmat, TMatrices, (-4,-3)),(-2,-1))/2
            elif opargs == 'σ_z' or opargs == 'σz':
                opmat = zflip_tyty(lMax)
                TMatrices = apply_ndmatrix_left(opmat, apply_ndmatrix_left(opmat, TMatrices, (-4,-3)),(-2,-1))/2
            elif mCN:
                rotN = int(mCN.group(2))
                power = int(mCN.group(1)) if mCN.group(1) else 1
                rotangle = 2*np.pi*power / rotN
                rot = WignerD_yy_fromvector(lMax, np.array([0,0,rotangle]))
                rotinv = WignerD_yy_fromvector(lMax, np.array([0,0,-rotangle]))
                TMatrices = apply_matrix_left(rot, apply_matrix_left(rotinv, TMatrices, -3), -1)
            else:
                raise ValueError('not implemented: ', opargs)
        else:
            raise ValueError('not implemented: ', optype)
    return (TMatrices, freqs, lMax)

class TMatrix(TMatrixSpec):
    '''
    TODO doc

    TODO support for different/multiple interpolators
    '''
    def __init__(self, tmatrix_spec):
        #self.specification = tmatrix_spec
        self.lMax_override = tmatrix_spec.lMax_override
        self.tmatrix_path = tmatrix_spec.tmatrix_path
        self.ops = tmatrix_spec.ops
        self.tmdata, self.freqs, self.lMax = get_TMatrix_fromspec(tmatrix_spec)
        self.nelem = get_nelem(self.lMax) 
        #self._interpolators = dict()
        self.default_interpolator = interpolate.interp1d(self.freqs, 
            self.tmdata, axis=0, kind='linear', fill_value='extrapolate')
        self.normalization = NormalizationT.TAYLOR # TODO others are not supported by the loading functions
    
    def atfreq(self, freq):
        freqarray = np.array(freq, copy=False)
        if freqarray.shape: # not just a scalar
            tm_interp = np.empty(freqarray.shape + (2, self.nelem, 2, self.nelem), dtype=np.complex_)
            for i in np.ndindex(freqarray.shape):
                tm_interp[i] = self.default_interpolator(freqarray[i])
            return tm_interp
        else: # scalar
            return self.default_interpolator(freq)

    __getitem__ = atfreq # might be changed later, use atfreq to be sure

def perform_tmspecs(tmspecs):
    """Takes a sequence of TMatrixSpec or TMatrix instances and returns
    a list of corresponding TMatrix instances"""
    return [(tmspec if hasattr(tmspec, "tmdata") else TMatrix(tmspec))
                for tmspec in tmspecs]
     
