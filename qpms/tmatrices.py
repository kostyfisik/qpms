import numpy as np
import quaternion, spherical_functions as sf # because of the Wigner matrices. These imports are SLOW.

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
        return 0
    elif (EM == b'M'):
        return 1
    else:
        return None

def loadScuffTMatrices(fileName):
    """
    TODO doc
    """
    μm = 1e-6
    table = np.genfromtxt(fileName, 
                  converters={1: _scuffTMatrixConvert_EM_01, 4: _scuffTMatrixConvert_EM_01},
                  dtype=[('freq', '<f8'), ('outc_type', '<i8'), ('outc_l', '<i8'), ('outc_m', '<i8'), 
                         ('inc_type', '<i8'), ('inc_l', '<i8'), ('inc_m', '<i8'), ('Treal', '<f8'), ('Timag', '<f8')]
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
    for inc_type in [0,1]:
        for outc_type in [0,1]:
            TMatrices[:,1-outc_type,:,1-inc_type,:] = TMatrices_tmp_real[:,:,outc_type,:,inc_type]+1j*TMatrices_tmp_imag[:,:,outc_type,:,inc_type]
    # IMPORTANT: now we are going from Reid's/Kristensson's/Jackson's/whoseever convention to Taylor's convention
    TMatrices[:,:,:,:,:] = TMatrices[:,:,:,:,:] * np.sqrt(ny*(ny+1))[ň,ň,ň,ň,:] / np.sqrt(ny*(ny+1))[ň,ň,:,ň,ň]
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
