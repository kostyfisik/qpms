import numpy as np
from qpms_c import *
ň = np.newaxis
import scipy
from scipy.constants import epsilon_0 as ε_0, c, pi as π, e, hbar as ℏ, mu_0 as μ_0
eV = e
from scipy.special import lpmn, lpmv, sph_jn, sph_yn, poch, gammaln
from scipy.misc import factorial
import math
import cmath

"""
'''
Try to import numba. Its pre-0.28.0 versions can not handle functions
containing utf8 identifiers, so we keep track about that.
'''
try:
    import numba
    use_jit = True
    if numba.__version__ >= '0.28.0':
        use_jit_utf8 = True
    else:
        use_jit_utf8 = False
except ImportError:
    use_jit = False
    use_jit_utf8 = False
'''
Accordingly, we define our own jit decorator that handles
different versions of numba or does nothing if numba is not
present. Note that functions that include unicode identifiers 
must be decorated with #@ujit
'''
#def dummywrap(f):
#    return f
def jit(f):
    if use_jit:
        return numba.jit(f)
    else:
        return f
def ujit(f):
    if use_jit_utf8:
        return numba.jit(f)
    else:
        return f
"""

# Coordinate transforms for arrays of "arbitrary" shape
#@ujit
def cart2sph(cart,axis=-1, allow2d=False):
    cart = np.array(cart, copy=False)
    if cart.shape[axis] == 3:
        [x, y, z] = np.split(cart,3,axis=axis)    
        r = np.linalg.norm(cart,axis=axis,keepdims=True)
        r_zero = np.logical_not(r)
        θ = np.arccos(z/(r+r_zero))
    elif cart.shape[axis] == 2 and allow2d:
        [x, y] = np.split(cart,2,axis=axis)
        r = np.linalg.norm(cart,axis=axis,keepdims=True)
        r_zero = np.logical_not(r)
        θ = np.broadcast_to(np.pi/2, x.shape)
    else:
        raise ValueError("The converted array has to have dimension 3 "
                         "(or 2 if allow2d==True)"
                         " along the given axis, not %d" % cart.shape[axis])

    φ = np.arctan2(y,x) # arctan2 handles zeroes correctly itself
    return np.concatenate((r,θ,φ),axis=axis)

#@ujit
def sph2cart(sph, axis=-1):
    sph = np.array(sph, copy=False)
    if (sph.shape[axis] != 3):
        raise ValueError("The converted array has to have dimension 3"
                         " along the given axis")
    [r,θ,φ] = np.split(sph,3,axis=axis)
    sinθ = np.sin(θ)
    x = r * sinθ * np.cos(φ)
    y = r * sinθ * np.sin(φ)
    z = r * np.cos(θ)
    return np.concatenate((x,y,z),axis=axis)

#@ujit
def sph_loccart2cart(loccart, sph, axis=-1):
    """
    Transformation of vector specified in local orthogonal coordinates 
    (tangential to spherical coordinates – basis r̂,θ̂,φ̂) to global cartesian
    coordinates (basis x̂,ŷ,ẑ)
    SLOW FOR SMALL ARRAYS
    
    Parameters
    ----------
    loccart: ... TODO
        the transformed vector in the local orthogonal coordinates
        
    sph: ... TODO
        the point (in spherical coordinates) at which the locally
        orthogonal basis is evaluated
        
    Returns
    -------
    output: ... TODO
        The coordinates of the vector in global cartesian coordinates
    """
    loccart = np.array(loccart, copy=False)
    sph = np.array(sph, copy=False)
    if (loccart.shape[axis] != 3):
        raise ValueError("The converted array has to have dimension 3"
                         " along the given axis")
    [r,θ,φ] = np.split(sph,3,axis=axis)
    sinθ = np.sin(θ)
    cosθ = np.cos(θ)
    sinφ = np.sin(φ)
    cosφ = np.cos(φ)
    
    #x = r * sinθ * cosφ
    #y = r * sinθ * sinφ
    #z = r * cosθ
    r̂x = sinθ * cosφ
    r̂y = sinθ * sinφ
    r̂z = cosθ
    θ̂x = cosθ * cosφ
    θ̂y = cosθ * sinφ
    θ̂z = -sinθ
    φ̂x = -sinφ
    φ̂y = cosφ
    φ̂z = np.zeros(φ̂y.shape)
    r̂ = np.concatenate((r̂x,r̂y,r̂z),axis=axis)
    θ̂ = np.concatenate((θ̂x,θ̂y,θ̂z),axis=axis)
    φ̂ = np.concatenate((φ̂x,φ̂y,φ̂z),axis=axis)
    [inr̂,inθ̂,inφ̂] = np.split(loccart,3,axis=axis)
    out=inr̂*r̂+inθ̂*θ̂+inφ̂*φ̂
    return out

#@ujit
def sph_loccart_basis(sph, sphaxis=-1, cartaxis=None):
    """
    Returns the local cartesian basis in terms of global cartesian basis.
    sphaxis refers to the original dimensions
    TODO doc
    """
    if(cartaxis is None):
        cartaxis = sph.ndim # default to last axis
    [r,θ,φ] = np.split(sph,3,axis=sphaxis)
    sinθ = np.sin(θ)
    cosθ = np.cos(θ)
    sinφ = np.sin(φ)
    cosφ = np.cos(φ)
    
    #x = r * sinθ * cosφ
    #y = r * sinθ * sinφ
    #z = r * cosθ
    r̂x = sinθ * cosφ
    r̂y = sinθ * sinφ
    r̂z = cosθ
    θ̂x = cosθ * cosφ
    θ̂y = cosθ * sinφ
    θ̂z = -sinθ
    φ̂x = -sinφ
    φ̂y = cosφ
    φ̂z = np.zeros(φ̂y.shape)
    #r̂ = np.concatenate((r̂x,r̂y,r̂z),axis=axis)
    #θ̂ = np.concatenate((θ̂x,θ̂y,θ̂z),axis=axis)
    #φ̂ = np.concatenate((φ̂x,φ̂y,φ̂z),axis=axis)
    x = np.expand_dims(np.concatenate((r̂x,θ̂x,φ̂x), axis=sphaxis),axis=cartaxis)
    y = np.expand_dims(np.concatenate((r̂y,θ̂y,φ̂y), axis=sphaxis),axis=cartaxis)
    z = np.expand_dims(np.concatenate((r̂z,θ̂z,φ̂z), axis=sphaxis),axis=cartaxis)
    out = np.concatenate((x,y,z),axis=cartaxis)
    return out


#@jit
def nelem2lMax(nelem):
    """
    Auxiliary inverse function to nelem(lMax) = (lMax + 2) * lMax. Returns 0 if
    it nelem does not come from positive integer lMax.
    """
    lMax = round(math.sqrt(1+nelem) - 1)
    if ((lMax < 1) or ((lMax + 2) * lMax != nelem)):
        return 0
    else:
        return lMax



#@jit
def lpy(nmax, z):
    """
    Associated legendre function and its derivatative at z in the 'y-indexing'.
    (Without Condon-Shortley phase AFAIK.)
    NOT THOROUGHLY TESTED
    
    Parameters
    ----------
    
    nmax: int
        The maximum order to which the Legendre functions will be evaluated..
        
    z: float
        The point at which the Legendre functions are evaluated.
        
    output: (P_y, dP_y) TODO
        y-indexed legendre polynomials and their derivatives
    
    """
    pmn_plus, dpmn_plus = lpmn(nmax, nmax, z)
    pmn_minus, dpmn_minus = lpmn(-nmax, nmax, z)
    nelem = nmax * nmax + 2*nmax
    P_y = np.empty((nelem), dtype=np.float_)
    dP_y = np.empty((nelem), dtype=np.float_)
    mn_p_y, mn_n_y = get_y_mn_unsigned(nmax)
    mn_plus_mask = (mn_p_y >= 0)
    mn_minus_mask = (mn_n_y >= 0)
    #print( mn_n_y[mn_minus_mask])
    P_y[mn_p_y[mn_plus_mask]] = pmn_plus[mn_plus_mask]
    P_y[mn_n_y[mn_minus_mask]] = pmn_minus[mn_minus_mask]
    dP_y[mn_p_y[mn_plus_mask]] = dpmn_plus[mn_plus_mask]
    dP_y[mn_n_y[mn_minus_mask]] = dpmn_minus[mn_minus_mask]
    return (P_y, dP_y)

def vswf_yr(pos_sph,lMax,J=1):
    """
    Normalized vector spherical wavefunctions $\widetilde{M}_{mn}^{j}$,
    $\widetilde{N}_{mn}^{j}$ as in [1, (2.40)].
    
    Parameters
    ----------
    
    pos_sph : np.array(dtype=float, shape=(someshape,3))
        The positions where the spherical vector waves are to be
        evaluated. The last axis corresponds to the individual
        points (r,θ,φ). The radial coordinate r is dimensionless,
        assuming that it has already been multiplied by the
        wavenumber.
    
    nmax : int
        The maximum order to which the VSWFs are evaluated.
        
    Returns
    -------
    
    output : np.array(dtype=complex, shape=(someshape,nmax*nmax + 2*nmax,3))
        Spherical vector wave functions evaluated at pos_sph,
        in the local basis (r̂,θ̂,φ̂). The last indices correspond
        to m, n (in the ordering given by mnindex()), and basis 
        vector index, respectively.
        
    [1] Jonathan M. Taylor. Optical Binding Phenomena: Observations and
    Mechanisms.
    """
    # TODO TRUEVECTORIZE
    pos_sph = np.array(pos_sph, copy = False)
    nelem = (lMax + 2) * lMax
    M_y = np.empty(pos_sph.shape[:-1] + (nelem,3), dtype = np.complex_)
    N_y = np.empty(pos_sph.shape[:-1] + (nelem,3), dtype = np.complex_)
    for i in np.ndindex(pos_sph.shape[:-1]):
        M_y[i], N_y[i] = vswf_yr1(pos_sph[i], lMax, J) # non-vectorised function
    return (M_y, N_y)

from scipy.special import sph_jn, sph_yn
#@jit
def _sph_zn_1(n,z):
    return sph_jn(n,z)
#@jit
def _sph_zn_2(n,z):
    return sph_yn(n,z)
#@jit
def _sph_zn_3(n,z):
    besj=sph_jn(n,z)
    besy=sph_yn(n,z)
    return (besj[0] + 1j*besy[0],besj[1] + 1j*besy[1])
#@jit
def _sph_zn_4(n,z):
    besj=sph_jn(n,z)
    besy=sph_yn(n,z)
    return (besj[0] - 1j*besy[0],besj[1] - 1j*besy[1])
_sph_zn = [_sph_zn_1,_sph_zn_2,_sph_zn_3,_sph_zn_4]

# computes bessel/hankel functions for orders from 0 up to n; drops
# the derivatives which are also included in scipy.special.sph_jn/yn
#@jit
def zJn(n, z, J=1):
    return _sph_zn[J-1](n=n,z=z)



# The following 4 funs have to be refactored, possibly merged

# FIXME: this can be expressed simply as:
# $$ -\frac{1}{2}\sqrt{\frac{2n+1}{4\pi}n\left(n+1\right)}(\delta_{m,1}+\delta_{m,-1}) $$
#@ujit
def π̃_zerolim(nmax): # seems OK
    """
    lim_{θ→ 0-} π̃(cos θ)
    """
    my, ny = get_mn_y(nmax)
    nelems = len(my)
    π̃_y = np.zeros((nelems))
    plus1mmask = (my == 1)
    minus1mmask = (my == -1)
    pluslim = -ny*(1+ny)/2
    minuslim = 0.5
    π̃_y[plus1mmask] = pluslim[plus1mmask]
    π̃_y[minus1mmask] = - minuslim
    prenorm =  np.sqrt((2*ny + 1)*factorial(ny-my)/(4*π*factorial(ny+my)))
    π̃_y = prenorm *     π̃_y
    return π̃_y

#@ujit
def π̃_pilim(nmax): # Taky OK, jen to možná není kompatibilní se vzorečky z mathematiky
    """
    lim_{θ→ π+} π̃(cos θ)
    """
    my, ny = get_mn_y(nmax)
    nelems = len(my)
    π̃_y = np.zeros((nelems))
    plus1mmask = (my == 1)
    minus1mmask = (my == -1)
    pluslim = (-1)**ny*ny*(1+ny)/2
    minuslim = 0.5*(-1)**ny
    π̃_y[plus1mmask] = pluslim[plus1mmask]
    π̃_y[minus1mmask] = minuslim[minus1mmask]
    prenorm =  np.sqrt((2*ny + 1)*factorial(ny-my)/(4*π*factorial(ny+my)))
    π̃_y = prenorm *     π̃_y
    return π̃_y

# FIXME: this can be expressed simply as
# $$ -\frac{1}{2}\sqrt{\frac{2n+1}{4\pi}n\left(n+1\right)}(\delta_{m,1}-\delta_{m,-1}) $$
#@ujit
def τ̃_zerolim(nmax):
    """
    lim_{θ→ 0-} τ̃(cos θ)
    """
    p0 = π̃_zerolim(nmax)
    my, ny = get_mn_y(nmax)
    minus1mmask = (my == -1)
    p0[minus1mmask] = -p0[minus1mmask]
    return p0

#@ujit
def τ̃_pilim(nmax):
    """
    lim_{θ→  π+} τ̃(cos θ)
    """
    t = π̃_pilim(nmax)
    my, ny = get_mn_y(nmax)
    plus1mmask = (my == 1)
    t[plus1mmask] = -t[plus1mmask]
    return t
    
#@ujit
def get_π̃τ̃_y1(θ,nmax):
    # TODO replace with the limit functions (below) when θ approaches
    # the extreme values at about 1e-6 distance
    """
    (... TODO)
    
    """
    if (abs(θ)<1e-6):
        return (π̃_zerolim(nmax),τ̃_zerolim(nmax))
    if (abs(θ-π)<1e-6):
        return (π̃_pilim(nmax),τ̃_pilim(nmax))
    my, ny = get_mn_y(nmax)
    nelems = len(my)
    Py, dPy = lpy(nmax, math.cos(θ))
    prenorm =  np.sqrt((2*ny + 1)*factorial(ny-my)/(4*π*factorial(ny+my)))
    π̃_y = prenorm * my * Py / math.sin(θ)  # bacha, možné dělení nulou
    τ̃_y = prenorm * dPy * (- math.sin(θ))  # TADY BACHA!!!!!!!!!! * (- math.sin(pos_sph[1])) ???
    return (π̃_y,τ̃_y)
    
#@ujit
def vswf_yr1(pos_sph,nmax,J=1):
    """
    As vswf_yr, but evaluated only at single position (i.e. pos_sph has
    to have shape=(3))
    """
    if (pos_sph[1].imag or pos_sph[2].imag):
        raise ValueError("The angles for the spherical wave functions can not be complex")
    kr = pos_sph[0] if pos_sph[0].imag else pos_sph[0].real # To supress the idiotic warning in scipy.special.sph_jn
    θ = pos_sph[1].real
    φ = pos_sph[2].real
    my, ny = get_mn_y(nmax)
    Py, dPy = lpy(nmax, math.cos(θ))
    nelems = nmax*nmax + 2*nmax
    # TODO Remove these two lines in production:
    if(len(Py) != nelems or len(my) != nelems):
        raise ValueError("This is very wrong.")
    prenorm =  np.sqrt((2*ny + 1)*factorial(ny-my)/(4*π*factorial(ny+my)))
    if (abs(θ)<1e-6): # Ošetření limitního chování derivací Leg. fcí
        π̃_y=π̃_zerolim(nmax)
        τ̃_y=τ̃_zerolim(nmax)
    elif (abs(θ-π)<1e-6):
        π̃_y=π̃_pilim(nmax)
        τ̃_y=τ̃_pilim(nmax)
    else:
        π̃_y = prenorm * my * Py / math.sin(θ) 
        τ̃_y = prenorm * dPy * (- math.sin(θ))  # TADY BACHA!!!!!!!!!! * (- math.sin(pos_sph[1])) ???
    z_n, dz_n = zJn(nmax, kr, J=J)
    z_y = z_n[ny]
    dz_y = dz_n[ny]
    eimf_y = np.exp(1j*my*φ) # zbytečné opakování my, lepší by bylo to spočítat jednou a vyindexovat
    M̃_y = np.zeros((nelems,3), dtype=np.complex_)
    M̃_y[:,1] = 1j * π̃_y * eimf_y * z_y
    M̃_y[:,2] =  -   τ̃_y * eimf_y * z_y
    Ñ_y = np.empty((nelems,3), dtype=np.complex_)
    Ñ_y[:,0] = (ny*(ny+1)/kr) * prenorm * Py * eimf_y * z_y
    Ñradial_fac_y = z_y / kr + dz_y
    Ñ_y[:,1] =    τ̃_y * eimf_y * Ñradial_fac_y
    Ñ_y[:,2] = 1j*π̃_y * eimf_y * Ñradial_fac_y
    return(M̃_y, Ñ_y)
    
#def plane_E_y(nmax):
#    """
#    The E_mn normalization factor as in [1, (3)] WITHOUT the E_0 factor,
#    y-indexed
#    
#    (... TODO)
#    
#    References
#    ----------
#    [1] Jonathan M. Taylor. Optical Binding Phenomena: Observations and
#    Mechanisms. FUCK, I MADE A MISTAKE: THIS IS FROM 7U
#    """
#    my, ny = get_mn_y(nmax)
#    return 1j**ny * np.sqrt((2*ny+1)*factorial(ny-my) /
#                            (ny*(ny+1)*factorial(ny+my))
#    )
#@ujit
def zplane_pq_y(nmax, betap = 0):
    """
    The z-propagating plane wave expansion coefficients as in [1, (1.12)]

    Taylor's normalization
    
    (... TODO)
    """
    my, ny = get_mn_y(nmax)
    U_y = 4*π * 1j**ny / (ny * (ny+1))
    π̃_y = π̃_zerolim(nmax)
    τ̃_y = τ̃_zerolim(nmax)
    
    # fixme co je zač ten e_θ ve vzorečku? (zde neimplementováno)
    p_y = U_y*(τ̃_y*math.cos(betap) - 1j*math.sin(betap)*π̃_y)
    q_y = U_y*(π̃_y*math.cos(betap) - 1j*math.sin(betap)*τ̃_y)
    return (p_y, q_y)
    
    
#import warnings
#@ujit
def plane_pq_y(nmax, kdir_cart, E_cart):
    """
    The plane wave expansion coefficients for any direction kdir_cart
    and amplitude vector E_cart (which might be complex, depending on
    the phase and polarisation state). If E_cart and kdir_cart are
    not orthogonal, the result should correspond to the k-normal part
    of E_cart.

    The Taylor's convention on the expansion is used; therefore,
    the electric field is expressed as
    E(r) = E_cart * exp(1j * k ∙ r)
         = -1j * ∑_y ( p_y * Ñ_y(|k| * r) + q_y * M̃_y(|k| * r))
    (note the -1j factor).

    Parameters
    ----------
    nmax: int
        The order of the expansion.
    kdir_cart: (...,3)-shaped real array
        The wave vector (its magnitude does not play a role).
    E_cart: (...,3)-shaped complex array
        Amplitude of the plane wave

    Returns
    -------
    p_y, q_y:
        The expansion coefficients for the electric (Ñ) and magnetic
        (M̃) waves, respectively.
    """
    # TODO TRUEVECTORIZE
    if np.iscomplexobj(kdir_cart):
        warnings.warn("The direction vector for the plane wave coefficients should be real. I am discarding the imaginary part now.")
        kdir_cart = kdir_cart.real
    
    E_cart = np.array(E_cart, copy=False)
    k_sph = cart2sph(kdir_cart)
    my, ny = get_mn_y(nmax)
    nelem = len(my)
    U_y = 4*π * 1j**ny / (ny * (ny+1))
    θ̂ = sph_loccart2cart(np.array([0,1,0]), k_sph, axis=-1)
    φ̂ = sph_loccart2cart(np.array([0,0,1]), k_sph, axis=-1)
    
    # not properly vectorised part:
    π̃_y = np.empty(k_sph.shape[:-1] + (nelem,), dtype=np.complex_)
    τ̃_y = np.empty(π̃_y.shape, dtype=np.complex_)
    for i in np.ndindex(k_sph.shape[:-1]):
        π̃_y[i], τ̃_y[i] = get_π̃τ̃_y1(k_sph[i][1], nmax) # this shit is not vectorised

    # last indices of the summands: [...,y,local cartesian component index (of E)]
    p_y = np.sum( U_y[...,:,ň]
                  * np.conj(np.exp(1j*my[...,:,ň]*k_sph[...,ň,ň,2]) * (
                     θ̂[...,ň,:]*τ̃_y[...,:,ň] + 1j*φ̂[...,ň,:]*π̃_y[...,:,ň]))
                  * E_cart[...,ň,:],
              axis=-1)
    q_y = np.sum( U_y[...,:,ň]
                  * np.conj(np.exp(1j*my[...,:,ň]*k_sph[...,ň,ň,2]) * (
                     θ̂[...,ň,:]*π̃_y[...,:,ň] + 1j*φ̂[...,ň,:]*τ̃_y[...,:,ň]))
                  * E_cart[...,ň,:],
          axis=-1)
    return (p_y, q_y)
    



# Material parameters
#@ujit
def ε_drude(ε_inf, ω_p, γ_p, ω): # RELATIVE permittivity, of course
    return ε_inf - ω_p*ω_p/(ω*(ω+1j*γ_p))


# In[8]:

# Mie scattering
#@ujit
def mie_coefficients(a, nmax,  #ω, ε_i, ε_e=1, J_ext=1, J_scat=3
                               k_i, k_e, μ_i=1, μ_e=1, J_ext=1, J_scat=3):
    """

    FIXME test the magnetic case
    TODO description
    RH concerns the N ("electric") part, RV the M ("magnetic") part
    #
    
    Parameters
    ----------
    a : float
        Diameter of the sphere.
        
    nmax : int
        To which order (inc. nmax) to compute the coefficients.
    
    ω : float
        Frequency of the radiation
    
    ε_i, ε_e, μ_i, μ_e : complex
        Relative permittivities and permeabilities of the sphere (_i)
        and the environment (_e)
        MAGNETIC (μ_i, μ_e != 1)  CASE UNTESTED AND PROBABLY BUGGY
    
    J_ext, J_scat : 1, 2, 3, or 4 (must be different)
        Specifies the species of the Bessel/Hankel functions in which
        the external incoming (J_ext) and scattered (J_scat) fields
        are represented. 1,2,3,4 correspond to j,y,h(1),h(2), respectively.
        The returned coefficients are always with respect to the decomposition
        of the "external incoming" wave.
    
    Returns
    -------
    RV == a/p, RH == b/q, TV = d/p, TH = c/q
    TODO 
    what does it return on index 0???
    FIXME permeabilities
    """
    # permittivities are relative!
    # cf. worknotes
    #print("a, nmax, ε_m, ε_b, ω",a, nmax, ε_m, ε_b, ω)
    #k_i = cmath.sqrt(ε_i*μ_i) * ω / c
    x_i = k_i * a
    #k_e = cmath.sqrt(ε_e*μ_e) * ω / c
    x_e = k_e * a
    #print("Mie: phase at radius: x_i",x_i,"x_e",x_e)
    m = k_i/k_e#cmath.sqrt(ε_i*μ_i/(ε_e*μ_e))
    # We "need" the absolute permeabilities for the final formula
    # This is not the absolute wave impedance, because only their ratio
    # ηi/ηe is important for getting the Mie coefficients.
    η_inv_i = k_i / μ_i
    η_inv_e = k_e / μ_e
    #print("k_m, x_m,k_b,x_b",k_m, x_m,k_b,x_b)
    zi, ži = zJn(nmax, x_i, J=1)
    #Pi = (zi * x_i)
    #Di = (zi + x_i * ži) / Pi # Vzoreček Taylor (2.9)
    #ži = zi + x_i * ži
    ze, že = zJn(nmax, x_e, J=J_ext)
    #Pe = (ze * x_e)
    #De = (ze + x_e * že) / Pe # Vzoreček Taylor (2.9)
    #že = ze + x_e * že
    zs, žs = zJn(nmax, x_e, J=J_scat)
    #Ps = (zs * x_e)
    #Ds = (zs + x_e * žs) / Ps # Vzoreček Taylor (2.9)
    #žs = zs + x_e * zs
    #RH = (μ_i*zi*že - μ_e*ze*ži) / (μ_i*zi*žs - μ_e*zs*ži)
    #RV = (μ_e*m*m*zi*že - μ_i*ze*ži) / (μ_e*m*m*zi*žs - μ_i*zs*ži)
    #TH = (μ_i*ze*žs - μ_i*zs*že) / (μ_i*zi*žs - μ_e*zs*ži)
    #TV = (μ_i*m*ze*žs - μ_i*m*zs*že) / (μ_e*m*m*zi*žs - μ_i*zs*ži)
    ži = zi/x_i+ži
    žs = zs/x_e+žs
    že = ze/x_e+že
    RV = -((-η_inv_i * že * zi + η_inv_e * ze * ži)/(-η_inv_e * ži * zs + η_inv_i * zi * žs))
    RH = -((-η_inv_e * že * zi + η_inv_i * ze * ži)/(-η_inv_i * ži * zs + η_inv_e * zi * žs))
    TV = -((-η_inv_e * že * zs + η_inv_e * ze * žs)/( η_inv_e * ži * zs - η_inv_i * zi * žs))
    TH = -(( η_inv_e * že * zs - η_inv_e * ze * žs)/(-η_inv_i * ži * zs + η_inv_e * zi * žs)) 
    return (RH, RV, TH, TV)

#@ujit
def G_Mie_scat_precalc_cart_new(source_cart, dest_cart, RH, RV, a, nmax, k_i, k_e, μ_i=1, μ_e=1, J_ext=1, J_scat=3):
    """
    Implementation according to Kristensson, page 50
    My (Taylor's) basis functions are normalized to n*(n+1), whereas Kristensson's to 1
    TODO: check possible -1 factors (cf. Kristensson's dagger notation)
    """
    my, ny = get_mn_y(nmax)
    nelem = len(my)
    #source to origin
    source_sph = cart2sph(source_cart)
    source_sph[0] = k_e * source_sph[0] 
    dest_sph = cart2sph(dest_cart)
    dest_sph[0] = k_e * dest_sph[0]
    if(dest_sph[0].real >= source_sph[0].real):
        lo_sph = source_sph
        hi_sph = dest_sph
    else:
        lo_sph = dest_sph
        hi_sph = source_sph
    lo_sph = source_sph
    hi_sph = dest_sph
        
    M̃lo_y, Ñlo_y = vswf_yr1(lo_sph,nmax,J=J_scat)
    lo_loccart_basis = sph_loccart_basis(lo_sph, sphaxis=-1, cartaxis=None)
    M̃lo_cart_y = np.sum(M̃lo_y[:,:,ň]*lo_loccart_basis[ň,:,:],axis=-2)
    Ñlo_cart_y = np.sum(Ñlo_y[:,:,ň]*lo_loccart_basis[ň,:,:],axis=-2)
    
    M̃hi_y, Ñhi_y = vswf_yr1(hi_sph,nmax,J=J_scat)#J_scat
    hi_loccart_basis = sph_loccart_basis(hi_sph, sphaxis=-1, cartaxis=None)
    M̃hi_cart_y = np.sum(M̃hi_y[:,:,ň]*hi_loccart_basis[ň,:,:],axis=-2)
    Ñhi_cart_y = np.sum(Ñhi_y[:,:,ň]*hi_loccart_basis[ň,:,:],axis=-2)
    
    G_y = (RH[ny][:,ň,ň] * M̃lo_cart_y[:,:,ň].conj() * M̃hi_cart_y[:,ň,:] + 
           RV[ny][:,ň,ň] * Ñlo_cart_y[:,:,ň].conj() * Ñhi_cart_y[:,ň,:]) / (ny * (ny+1))[:,ň,ň]
    return 1j* k_e*np.sum(G_y,axis=0)
    

# From PRL 112, 253601 (1)
#@ujit
def Grr_Delga(nmax, a, r, k, ε_m, ε_b):
    om = k * c
    z = (r-a)/a
    g0 = om*cmath.sqrt(ε_b)/(6*c*π)
    n = np.arange(1,nmax+1)
    s = np.sum( (n+1)**2 * (ε_m-ε_b) / ((1+z)**(2*n+4) * (ε_m + ((n+1)/n)*ε_b)))
    return (g0 + s * c**2/(4*π*om**2*ε_b*a**3))
    


# TODOs
# ====
# 
# Rewrite the functions zJn, lpy in (at least simulated) universal manner.
# Then universalise the rest
# 
# Implement the actual multiple scattering
# 
# Test if the decomposition of plane wave works also for absorbing environment (complex k).

# From PRL 112, 253601 (1)
#@ujit
def Grr_Delga(nmax, a, r, k, ε_m, ε_b):
    om = k * c
    z = (r-a)/a
    g0 = om*cmath.sqrt(ε_b)/(6*c*π)
    n = np.arange(1,nmax+1)
    s = np.sum( (n+1)**2 * (ε_m-ε_b) / ((1+z)**(2*n+4) * (ε_m + ((n+1)/n)*ε_b)))
    return (g0 + s * c**2/(4*π*om**2*ε_b*a**3))

#@ujit
def G0_dip_1(r_cart,k):
    """
    Free-space dyadic Green's function in terms of the spherical vector waves.
    FIXME
    """
    sph = cart2sph(r_cart*k)
    pfz = 0.32573500793527994772 # 1./math.sqrt(3.*π)
    pf = 0.23032943298089031951 # 1./math.sqrt(6.*π)
    M1_y, N1_y = vswf_yr1(sph,nmax = 1,J=3)
    loccart_basis = sph_loccart_basis(sph, sphaxis=-1, cartaxis=None)
    N1_cart = np.sum(N1_y[:,:,ň]*loccart_basis[ň,:,:],axis=-2)
    coeffs_cart = np.array([[pf,-1j*pf,0.],[0.,0.,pfz],[-pf,-1j*pf,0.]]).conj()
    return 1j*k*np.sum(coeffs_cart[:,:,ň]*N1_cart[:,ň,:],axis=0)/2.

# Free-space dyadic Green's functions from RMP 70, 2, 447 =: [1]
# (The numerical value is correct only at the regular part, i.e. r != 0)
#@ujit
def _P(z):
    return (1-1/z+1/(z*z))
#@ujit
def _Q(z):
    return (-1+3/z-3/(z*z))

# [1, (9)] FIXME The sign here is most likely wrong!!!
#@ujit
def G0_analytical(r #cartesian!
                  , k):
    I=np.identity(3)
    rn = sph_loccart2cart(np.array([1.,0.,0.]), cart2sph(r), axis=-1)
    rnxrn = rn[...,:,ň] * rn[...,ň,:]
    r = np.linalg.norm(r, axis=-1)
    #print(_P(1j*k*r).shape,_Q(1j*k*r).shape, rnxrn.shape, I.shape)
    return ((-np.exp(1j*k*r)/(4*π*r))[...,ň,ň] *
                   (_P(1j*k*r)[...,ň,ň]*I
                    +_Q(1j*k*r)[...,ň,ň]*rnxrn
                   ))

# [1, (11)]
#@ujit
def G0L_analytical(r, k):
    I=np.identity(3)
    rn = sph_loccart2cart(np.array([1.,0.,0.]), cart2sph(r), axis=-1)
    rnxrn = rn[...,:,ň] * rn[...,ň,:]
    r = np.linalg.norm(r, axis=-1)
    return (I-3*rnxrn)/(4*π*k*k*r**3)[...,ň,ň]

# [1,(10)]
#@jit
def G0T_analytical(r, k):
    return G0_analytical(r,k) - G0L_analytical(r,k)

#@ujit
def G0_sum_1_slow(source_cart, dest_cart, k, nmax):
    my, ny = get_mn_y(nmax)
    nelem = len(my) 
    RH = np.full((nelem),1)
    RV = RH
    return G_Mie_scat_precalc_cart(source_cart, dest_cart, RH, RV, a=0.001, nmax=nmax, k_i=1, k_e=k, μ_i=1, μ_e=1, J_ext=1, J_scat=3)


