import numpy as np
from qpms_c import *
ň = np.newaxis
from scipy.constants import epsilon_0 as ε_0, c, pi as π, e, hbar as ℏ, mu_0 as μ_0
eV = e
from scipy.special import lpmn, lpmv, sph_jn, sph_yn, poch
from scipy.misc import factorial
import math
import cmath
import quaternion, spherical_functions as sf # because of the Wigner matrices


# Coordinate transforms for arrays of "arbitrary" shape
def cart2sph(cart,axis=-1):
    if (cart.shape[axis] != 3):
        raise ValueError("The converted array has to have dimension 3"
                         " along the given axis")
    [x, y, z] = np.split(cart,3,axis=axis)
    r = np.linalg.norm(cart,axis=axis,keepdims=True)
    r_zero = np.logical_not(r)
    θ = np.arccos(z/(r+r_zero))
    φ = np.arctan2(y,x) # arctan2 handles zeroes correctly itself
    return np.concatenate((r,θ,φ),axis=axis)

def sph2cart(sph, axis=-1):
    if (sph.shape[axis] != 3):
        raise ValueError("The converted array has to have dimension 3"
                         " along the given axis")
    [r,θ,φ] = np.split(sph,3,axis=axis)
    sinθ = np.sin(θ)
    x = r * sinθ * np.cos(φ)
    y = r * sinθ * np.sin(φ)
    z = r * np.cos(θ)
    return np.concatenate((x,y,z),axis=axis)

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

def vswf_yr(pos_sph,nmax,J=1):
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
    #mi, ni = mnindex(nmax)
    #nelems = nmax*nmax + 2*nmax
    ## TODO Remove these two lines in production:
    #if(len(mi) != nelems):
    #    raise ValueError("This is very wrong.")
    ## Pre-calculate the associated Legendre function
    #Prmn, dPrmn = lpmn(nmax,nmax,)
    ## Normalized funs π̃, τ̃
    #π̃ = 
    pass

from scipy.special import sph_jn, sph_yn
def _sph_zn_1(n,z):
    return sph_jn(n,z)
def _sph_zn_2(n,z):
    return sph_yn(n,z)
def _sph_zn_3(n,z):
    besj=sph_jn(n,z)
    besy=sph_yn(n,z)
    return (besj[0] + 1j*besy[0],besj[1] + 1j*besy[1])
def _sph_zn_4(n,z):
    besj=sph_jn(n,z)
    besy=sph_yn(n,z)
    return (besj[0] - 1j*besy[0],besj[1] - 1j*besy[1])
_sph_zn = [_sph_zn_1,_sph_zn_2,_sph_zn_3,_sph_zn_4]

# computes bessel/hankel functions for orders from 0 up to n; drops
# the derivatives which are also included in scipy.special.sph_jn/yn
def zJn(n, z, J=1):
    return _sph_zn[J-1](n=n,z=z)



# The following 4 funs have to be refactored, possibly merged

# FIXME: this can be expressed simply as:
# $$ -\frac{1}{2}\sqrt{\frac{2n+1}{4\pi}n\left(n+1\right)}(\delta_{m,1}+\delta_{m,-1}) $$
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
def τ̃_zerolim(nmax):
    """
    lim_{θ→ 0-} τ̃(cos θ)
    """
    p0 = π̃_zerolim(nmax)
    my, ny = get_mn_y(nmax)
    minus1mmask = (my == -1)
    p0[minus1mmask] = -p0[minus1mmask]
    return p0

def τ̃_pilim(nmax):
    """
    lim_{θ→  π+} τ̃(cos θ)
    """
    t = π̃_pilim(nmax)
    my, ny = get_mn_y(nmax)
    plus1mmask = (my == 1)
    t[plus1mmask] = -t[plus1mmask]
    return t
    
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

def zplane_pq_y(nmax, betap = 0):
    """
    The z-propagating plane wave expansion coefficients as in [1, (1.12)]
    
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
def plane_pq_y(nmax, kdir_cart, E_cart):
    """
    The plane wave expansion coefficients for any direction kdir_cart
    and amplitude vector E_cart (which might be complex, depending on
    the phase and polarisation state). If E_cart and kdir_cart are
    not orthogonal, the result should correspond to the k-normal part
    of E_cart.
    """
    if np.iscomplexobj(kdir_cart):
        warnings.warn("The direction vector for the plane wave coefficients should be real. I am discarding the imaginary part now.")
        kdir_cart = kdir_cart.real
        
    k_sph = cart2sph(kdir_cart)
    π̃_y, τ̃_y = get_π̃τ̃_y1(k_sph[1], nmax) 
    my, ny = get_mn_y(nmax)
    U_y = 4*π * 1j**ny / (ny * (ny+1))
    θ̂ = sph_loccart2cart(np.array([0,1,0]), k_sph, axis=-1)
    φ̂ = sph_loccart2cart(np.array([0,0,1]), k_sph, axis=-1)
    p_y = np.sum( U_y[:,ň]
                  * np.conj(np.exp(1j*my[:,ň]*k_sph[2]) * (
                     θ̂[ň,:]*τ̃_y[:,ň] + 1j*φ̂[ň,:]*π̃_y[:,ň]))
                  * E_cart[ň,:],
              axis=-1)
    q_y = np.sum( U_y[:,ň]
                  * np.conj(np.exp(1j*my[:,ň]*k_sph[2]) * (
                     θ̂[ň,:]*π̃_y[:,ň] + 1j*φ̂[ň,:]*τ̃_y[:,ň]))
                  * E_cart[ň,:],
          axis=-1)
    return (p_y, q_y)
    


# Functions copied from scattering_xu, additionaly normalized
from py_gmm.gmm import vec_trans as vc

def q_max(m,n,μ,ν):
    return min(n,ν,(n+ν-abs(m+μ))/2)
    
# returns array with indices corresponding to q
# argument q does nothing for now
def a_q(m,n,μ,ν,q = None):
    qm=q_max(m,n,μ,ν)
    res, err= vc.gaunt_xu(m,n,μ,ν,qm)
    if(err):
        print("m,n,μ,ν,qm = ",m,n,μ,ν,qm)
        raise ValueError('Something bad in the fortran subroutine gaunt_xu happened')
    return res

# All arguments are single numbers (for now)
# ZDE VYCHÁZEJÍ DIVNÁ ZNAMÉNKA
def Ã(m,n,μ,ν,kdlj,θlj,φlj,r_ge_d,J):
    exponent=(math.lgamma(2*n+1)-math.lgamma(n+2)+math.lgamma(2*ν+3)-math.lgamma(ν+2) 
                +math.lgamma(n+ν+m-μ+1)-math.lgamma(n-m+1)-math.lgamma(ν+μ+1)
                +math.lgamma(n+ν+1) - math.lgamma(2*(n+ν)+1))
    presum = math.exp(exponent)
    presum = presum * np.exp(1j*(μ-m)*φlj) * (-1)**m * 1j**(ν+n) / (4*n)
    qmax = math.floor(q_max(-m,n,μ,ν)) #nemá tu být +m?
    q = np.arange(qmax+1, dtype=int)
    # N.B. -m !!!!!!
    a1q = a_q(-m,n,μ,ν) # there is redundant calc. of qmax
    ã1q = a1q / a1q[0]
    p = n+ν-2*q
    if(r_ge_d):
        J = 1
    zp = zJn(n+ν,kdlj,J)[0][p]
    Pp = lpmv(μ-m,p,math.cos(θlj))
    summandq = (n*(n+1) + ν*(ν+1) - p*(p+1)) * (-1)**q * ã1q * zp * Pp
  
    # Taylor normalisation v2, proven to be equivalent (NS which is better)
    prenormratio = 1j**(ν-n) * math.sqrt(((2*ν+1)/(2*n+1))* math.exp(
        math.lgamma(n+m+1)-math.lgamma(n-m+1)+math.lgamma(ν-μ+1)-math.lgamma(ν+μ+1)))
    presum = presum / prenormratio
    
    # Taylor normalisation
    #prenormmn =  math.sqrt((2*n + 1)*math.factorial(n-m)/(4*π*factorial(n+m)))
    #prenormμν =  math.sqrt((2*ν + 1)*math.factorial(ν-μ)/(4*π*factorial(ν+μ)))
    #presum = presum * prenormμν / prenormmn
    
    return presum * np.sum(summandq)
    
# ZDE OPĚT JINAK ZNAMÉNKA než v Xu (J. comp. phys 127, 285)
def B̃(m,n,μ,ν,kdlj,θlj,φlj,r_ge_d,J):
    exponent=(math.lgamma(2*n+3)-math.lgamma(n+2)+math.lgamma(2*ν+3)-math.lgamma(ν+2) 
                +math.lgamma(n+ν+m-μ+2)-math.lgamma(n-m+1)-math.lgamma(ν+μ+1)
                +math.lgamma(n+ν+2) - math.lgamma(2*(n+ν)+3))
    presum = math.exp(exponent)
    presum = presum * np.exp(1j*(μ-m)*φlj) * (-1)**m * 1j**(ν+n+1) / (
        (4*n)*(n+1)*(n+m+1))
    Qmax = math.floor(q_max(-m,n+1,μ,ν))
    q = np.arange(Qmax+1, dtype=int)
    if (μ == ν): # it would disappear in the sum because of the factor (ν-μ) anyway
        ã2q = 0
    else:
        a2q = a_q(-m-1,n+1,μ+1,ν)
        ã2q = a2q / a2q[0]
    a3q = a_q(-m,n+1,μ,ν)
    ã3q = a3q / a3q[0]
    #print(len(a2q),len(a3q))
    p = n+ν-2*q
    if(r_ge_d):
        J = 1
    zp_ = zJn(n+1+ν,kdlj,J)[0][p+1] # je ta +1 správně?
    Pp_ = lpmv(μ-m,p+1,math.cos(θlj))
    summandq = ((2*(n+1)*(ν-μ)*ã2q
                 -(-ν*(ν+1) - n*(n+3) - 2*μ*(n+1)+p*(p+3))* ã3q)
                *(-1)**q * zp_ * Pp_)
    
    # Taylor normalisation v2, proven to be equivalent
    prenormratio = 1j**(ν-n) * math.sqrt(((2*ν+1)/(2*n+1))* math.exp(
        math.lgamma(n+m+1)-math.lgamma(n-m+1)+math.lgamma(ν-μ+1)-math.lgamma(ν+μ+1)))
    presum = presum / prenormratio
    
    ## Taylor normalisation
    #prenormmn =  math.sqrt((2*n + 1)*math.factorial(n-m)/(4*π*factorial(n+m)))
    #prenormμν =  math.sqrt((2*ν + 1)*math.factorial(ν-μ)/(4*π*factorial(ν+μ)))
    #presum = presum * prenormμν / prenormmn
    
    return presum * np.sum(summandq)
 


# In[7]:

# Material parameters
def ε_drude(ε_inf, ω_p, γ_p, ω): # RELATIVE permittivity, of course
    return ε_inf - ω_p*ω_p/(ω*(ω+1j*γ_p))


# In[8]:

# Mie scattering
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
    
def G_Mie_scat_precalc_cart(source_cart, dest_cart, RH, RV, a, nmax, k_i, k_e, μ_i=1, μ_e=1, J_ext=1, J_scat=3):
    """
    r1_cart (destination), r2_cart (source) and the result are in cartesian coordinates
    the result indices are in the source-destination order
    TODO
    """
    my, ny = get_mn_y(nmax)
    nelem = len(my)
    #source to origin
    so_sph = cart2sph(-source_cart)
    kd_so = k_e * so_sph[0]
    θ_so = so_sph[1]
    φ_so = so_sph[2]
    # Decomposition of the source N_0,1, N_-1,1, and N_1,1 in the nanoparticle center
    p_0 = np.empty((nelem), dtype=np.complex_)
    q_0 = np.empty((nelem), dtype=np.complex_)
    p_minus = np.empty((nelem), dtype=np.complex_)
    q_minus = np.empty((nelem), dtype=np.complex_)
    p_plus = np.empty((nelem), dtype=np.complex_)
    q_plus = np.empty((nelem), dtype=np.complex_)
    for y in range(nelem):
        m = my[y]
        n = ny[y]
        p_0[y]     = Ã(m,n, 0,1,kd_so,θ_so,φ_so,False,J=J_scat)
        q_0[y]     = B̃(m,n, 0,1,kd_so,θ_so,φ_so,False,J=J_scat)
        p_minus[y] = Ã(m,n,-1,1,kd_so,θ_so,φ_so,False,J=J_scat)
        q_minus[y] = B̃(m,n,-1,1,kd_so,θ_so,φ_so,False,J=J_scat)
        p_plus[y]  = Ã(m,n, 1,1,kd_so,θ_so,φ_so,False,J=J_scat)
        q_plus[y]  = B̃(m,n, 1,1,kd_so,θ_so,φ_so,False,J=J_scat)
    a_0 = RV[ny] * p_0
    b_0 = RH[ny] * q_0
    a_plus = RV[ny] * p_plus
    b_plus = RH[ny] * q_plus
    a_minus = RV[ny] * p_minus
    b_minus = RH[ny] * q_minus
    orig2dest_sph = cart2sph(dest_cart)
    orig2dest_sph[0] = k_e*orig2dest_sph[0]
    M_dest_y, N_dest_y = vswf_yr1(orig2dest_sph,nmax,J=J_scat)
    # N.B. these are in the local cartesian coordinates (r̂,θ̂,φ̂)
    N_dest_0     = np.sum(a_0[:,ň]    * N_dest_y, axis=-2)
    M_dest_0     = np.sum(b_0[:,ň]    * M_dest_y, axis=-2)
    N_dest_plus  = np.sum(a_plus[:,ň] * N_dest_y, axis=-2)
    M_dest_plus  = np.sum(b_plus[:,ň] * M_dest_y, axis=-2)
    N_dest_minus = np.sum(a_minus[:,ň]* N_dest_y, axis=-2)
    M_dest_minus = np.sum(b_minus[:,ň]* M_dest_y, axis=-2)
    prefac = math.sqrt(1/(4*3*π))#/ε_0
    G_sourcez_dest = prefac * (N_dest_0+M_dest_0)
    G_sourcex_dest = prefac * (N_dest_minus+M_dest_minus-N_dest_plus-M_dest_plus)/math.sqrt(2)
    G_sourcey_dest = prefac * (N_dest_minus+M_dest_minus+N_dest_plus+M_dest_plus)/(1j*math.sqrt(2))
    G_source_dest = np.array([G_sourcex_dest, G_sourcey_dest, G_sourcez_dest])
    # To global cartesian coordinates:
    G_source_dest = sph_loccart2cart(G_source_dest, sph=orig2dest_sph, axis=-1)
    return G_source_dest
    
def G_Mie_scat_cart(source_cart, dest_cart, a, nmax, k_i, k_e, μ_i=1, μ_e=1, J_ext=1, J_scat=3):
    """
    TODO
    """
    RH, RV, TH, TV = mie_coefficients(a=a, nmax=nmax, k_i=k_i, k_e=k_e, μ_i=μ_i, μ_e=μ_e, J_ext=J_ext, J_scat=J_scat)
    return G_Mie_scat_precalc_cart_new(source_cart, dest_cart, RH, RV, a, nmax, k_i, k_e, μ_i, μ_e, J_ext, J_scat)


#TODO
def cross_section_Mie_precalc():
    pass

def cross_section_Mie(a, nmax, k_i, k_e, μ_i, μ_e,):
    pass


# In[9]:

# From PRL 112, 253601 (1)
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
def Grr_Delga(nmax, a, r, k, ε_m, ε_b):
    om = k * c
    z = (r-a)/a
    g0 = om*cmath.sqrt(ε_b)/(6*c*π)
    n = np.arange(1,nmax+1)
    s = np.sum( (n+1)**2 * (ε_m-ε_b) / ((1+z)**(2*n+4) * (ε_m + ((n+1)/n)*ε_b)))
    return (g0 + s * c**2/(4*π*om**2*ε_b*a**3))


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
def _P(z):
    return (1-1/z+1/(z*z))
def _Q(z):
    return (-1+3/z-3/(z*z))

# [1, (9)] FIXME The sign here is most likely wrong!!!
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
def G0L_analytical(r, k):
    I=np.identity(3)
    rn = sph_loccart2cart(np.array([1.,0.,0.]), cart2sph(r), axis=-1)
    rnxrn = rn[...,:,ň] * rn[...,ň,:]
    r = np.linalg.norm(r, axis=-1)
    return (I-3*rnxrn)/(4*π*k*k*r**3)[...,ň,ň]

# [1,(10)]
def G0T_analytical(r, k):
    return G0_analytical(r,k) - G0L_analytical(r,k)


def G0_sum_1_slow(source_cart, dest_cart, k, nmax):
    my, ny = get_mn_y(nmax)
    nelem = len(my) 
    RH = np.full((nelem),1)
    RV = RH
    return G_Mie_scat_precalc_cart(source_cart, dest_cart, RH, RV, a=0.001, nmax=nmax, k_i=1, k_e=k, μ_i=1, μ_e=1, J_ext=1, J_scat=3)



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
    
def yflip_yy(lmax):
    """
    TODO doc
    """
    my, ny = get_mn_y(lmax)
    elems = xflip_yy(lmax)
    elems[(my % 2)==1] = elems[(my % 2)==1] * -1 # Obvious sign of tiredness (this is correct but ugly; FIXME)
    return elems

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
    TMatrices = np.empty((len(freqs),2,nelem,2,nelem),dtype=complex)
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


####################
# Array simulations
####################


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

def scatter_plane_wave(omega, epsilon_b, positions, Tmatrices, k_dirs, E_0s, #saveto = None
                      ):
    """
    Solves the plane wave linear scattering problem for a structure of "non-touching" particles
    for one frequency and arbitrary number K of incoming plane waves.
    
    Parameters
    ----------
    omega : positive number
        The frequency of the field.
    epsilon_b : complex number
        Permittivity of the background medium (which has to be isotropic).
    positions : (N,3)-shaped real array
        Cartesian positions of the particles.
    TMatrices : (N,2,nelem,2,nelem) or compatible
        The T-matrices in the "Taylor convention" describing the scattering on a single nanoparticle.
        If all the particles are identical and equally oriented, only one T-matrix can be given.
        nelems = (lMax + 2) * lMax, where lMax is the highest multipole order to which the scattering
        is calculated.
    k_dirs : (K,3)-shaped real array or compatible
        The direction of the incident field wave vector, normalized to one.
    E_0s : (K,3)-shaped complex array or compatible
        The electric intensity amplitude of the incident field.
        
    Returns
    -------
    ab : (K, N, 2, nelem)-shaped complex array
        The a (electric wave), b (magnetic wave) coefficients of the outgoing field for each particle
    # Fuck this, it will be wiser to make separate function to calculate those from ab:
    # sigma_xxx : TODO (K, 2, nelem)
    #    TODO partial (TODO which?) cross-section for each type of outgoing waves, summed over all
    #    nanoparticles (total cross section is given by the sum of this.)
    """
    nelem = TMatrices.shape[-1]
    if ((nelem != TMatrices.shape[-3]) or (2 != TMatrices.shape[-2]) or (2 != TMatrices.shape[-4])):
        raise ValueError('The T-matrices must be of shape (N, 2, nelem, 2, nelem) but are of shape %s' % (str(TMatrices.shape),))
    lMax = nelem2lMax(nelem)
    if not lMax:
        raise ValueError('The "nelem" dimension of T-matrix has invalid value (%d).' % nelem)
    # TODO perhaps more checks.
    raise Error('Not implemented.')
    pass

import warnings
def scatter_plane_wave_rectarray(omega, epsilon_b, xN, yN, xd, yd, TMatrices, k_dirs, E_0s, 
        return_pq_0 = False, return_pq= False, return_xy = False):
    """
    Solves the plane wave linear scattering problem for a rectangular array of particles
    for one frequency and arbitrary number K of incoming plane waves.
    
    Parameters
    ----------
    omega : positive number
        The frequency of the field.
    epsilon_b : complex number
        Permittivity of the background medium (which has to be isotropic).
    xN, yN : positive integers
        Particle numbers in the x and y dimensions
    xd, yd : positive numbers
        Periodicities in the x and y direction
    TMatrices : (xN, yN,2,nelem,2,nelem) or compatible or (2,nelem,2,nelem)
        The T-matrices in the "Taylor convention" describing the scattering on a single nanoparticle.
        If all the particles are identical and equally oriented, only one T-matrix can be given.
        nelems = (lMax + 2) * lMax, where lMax is the highest multipole order to which the scattering
        is calculated.
        Electric wave index is 0, magnetic wave index is 1.
    k_dirs : (K,3)-shaped real array or compatible
        The direction of the incident field wave vector, normalized to one.
    E_0s : (K,3)-shaped complex array or compatible
        The electric intensity amplitude of the incident field.
    return_pq_0 : bool
        Return also the multipole decomposition coefficients of the incoming plane wave.
    return_pq : bool NOT IMPLEMENTED
        Return also the multipole decomposition coefficients of the field incoming to each
        particle (inc. the field scattered from other particles.
    return_xy : bool
        Return also the cartesian x, y positions of the particles. 
        
    Returns
    -------
    ab : (K, xN, yN, 2, nelem)-shaped complex array
        The a (electric wave), b (magnetic wave) coefficients of the outgoing field for each particle.
        If none of return_pq or return_xy is set, the array is not enclosed in a tuple.
    pq_0 : (K, xN, yn, 2, nelem)-shaped complex array
        The p_0 (electric wave), b_0 (magnetic wave) coefficients of the incoming plane wave for each particle.
    pq : (K, xN, yN, 2, nelem)-shaped complex array NOT IMPLEMENTED
        The p (electric wave), q (magnetic wave) coefficients of the total exciting field
        for each particle (including the field scattered from other particles)
    x, y : (xN, yN)-shaped real array
        The x,y positions of the nanoparticles.
    """
    nelem = TMatrices.shape[-1]
    if ((nelem != TMatrices.shape[-3]) or (2 != TMatrices.shape[-2]) or (2 != TMatrices.shape[-4])):
        raise ValueError('The T-matrices must be of shape (N, 2, nelem, 2, nelem) but are of shape %s' % (str(TMatrices.shape),))
    lMax = nelem2lMax(nelem)
    if not lMax:
        raise ValueError('The "nelem" dimension of T-matrix has invalid value (%d).' % nelem)
    # TODO perhaps more checks.
    k_out = omega * math.sqrt(epsilon_b) / c # wave number
    my, ny = get_mn_y(lMax)
    N = yN * xN
    
    J_scat=3
    J_ext=1
    
    # Do something with this ugly indexing crap
    xind, yind = np.meshgrid(np.arange(xN),np.arange(yN), indexing='ij')
    xind = xind.flatten()
    yind = yind.flatten()
    xyind = np.stack((xind, yind, np.zeros((xind.shape),dtype=int)),axis=-1)
    cart_lattice=xyind * np.array([xd, yd, 0])
    x=cart_lattice[:,0]
    y=cart_lattice[:,1]
    xyind = xyind[:,0:2]
    
    # Lattice speedup
    Agrid = np.zeros((nelem, 2*xN, 2*yN, nelem),dtype=np.complex_)
    Bgrid = np.zeros((nelem, 2*xN, 2*yN, nelem),dtype=np.complex_)
    for yl in range(nelem): # source
        for xij in range(2*xN):
            for yij in range(2*yN):
                for yj in range(nelem): #dest
                  if((yij != yN) or (xij != xN)):
                    d_l2j = cart2sph(np.array([(xij-xN)*xd, (yij-yN)*yd, 0]))
                    Agrid[yj, xij, yij, yl] = Ã(my[yj],ny[yj],my[yl],ny[yl],kdlj=d_l2j[0]*k_out,θlj=d_l2j[1],φlj=d_l2j[2],r_ge_d=False,J=J_scat)
                    Bgrid[yj, xij, yij, yl] = B̃(my[yj],ny[yj],my[yl],ny[yl],kdlj=d_l2j[0]*k_out,θlj=d_l2j[1],φlj=d_l2j[2],r_ge_d=False,J=J_scat)

    # Translation coefficient matrix T
    transmat = np.zeros((xN* yN, 2, nelem, xN* yN, 2, nelem),dtype=np.complex_)
    for l in range(N):
        xil, yil = xyind[l]
        for j in range(N):
            xij, yij = xyind[j]
            if (l!=j):
                        transmat[j,0,:,l,0,:] = Agrid[:, xij - xil + xN, yij - yil + yN, :]
                        transmat[j,0,:,l,1,:] = Bgrid[:, xij - xil + xN, yij - yil + yN, :]
                        transmat[j,1,:,l,0,:] = Bgrid[:, xij - xil + xN, yij - yil + yN, :]
                        transmat[j,1,:,l,1,:] = Agrid[:, xij - xil + xN, yij - yil + yN, :]
    Agrid = None
    Bgrid = None

    # Now we solve a linear problem (1 - M T) A = M P_0 where M is the T-matrix :-)
    MT = np.empty((N,2,nelem,N,2,nelem),dtype=np.complex_)
    
    TMatrices = np.broadcast_to(TMatrices, (xN, yN, 2, nelem, 2, nelem))
    for j in range(N): # I wonder how this can be done without this loop...
        xij, yij = xyind[j]
        MT[j] = np.tensordot(TMatrices[xij, yij],transmat[j],axes=([-2,-1],[0,1]))
    MT.shape = (N*2*nelem, N*2*nelem)
    leftmatrix = np.identity(N*2*nelem) - MT
    MT = None

    if ((1 == k_dirs.ndim) and (1 == E_0s.ndim)):
        k_cart = k_dirs * k_out # wave vector of the incident plane wave
        pq_0 = np.zeros((N,2,nelem), dtype=np.complex_)
        p_y0, q_y0 = plane_pq_y(lMax, k_cart, E_0s)
        pq_0[:,0,:] = np.exp(1j*np.sum(k_cart[ň,:]*cart_lattice,axis=-1))[:, ň] * p_y0[ň, :]
        pq_0[:,1,:] = np.exp(1j*np.sum(k_cart[ň,:]*cart_lattice,axis=-1))[:, ň] * q_y0[ň, :]
        if (return_pq_0):
            pq_0_arr = pq_0
        MP_0 = np.empty((N,2,nelem),dtype=np.complex_)
        for j in range(N): # I wonder how this can be done without this loop...
            MP_0[j] = np.tensordot(TMatrices[xij, yij],pq_0[j],axes=([-2,-1],[-2,-1]))
        MP_0.shape = (N*2*nelem,)
        
        ab = np.linalg.solve(leftmatrix, MP_0)
        ab.shape = (xN, yN, 2, nelem)
    else:
        # handle "broadcasting" for k, E
        if 1 == k_dirs.ndim:
            k_dirs = k_dirs[ň,:]
        if 1 == E_0s.ndim:
            E_0s = E_0s[ň,:]
        K = max(E_0s.shape[-2], k_dirs.shape[-2])
        k_dirs = np.broadcast_to(k_dirs,(K,3))
        E_0s = np.broadcast_to(E_0s, (K,3))
        
        # А ну, чики-брики и в дамки!
        lupiv = scipy.linalg.lu_factor(leftmatrix, overwrite_a=True)
        leftmatrix = None
        
        if (return_pq_0):
            pq_0_arr = np.zeros((K,N,2,nelem), dtype=np.complex_)
        ab = np.empty((K,N*2*nelem), dtype=complex)
        for ki in range(K):
            k_cart = k_dirs[ki] * k_out
            pq_0 = np.zeros((N,2,nelem), dtype=np.complex_)
            p_y0, q_y0 = plane_pq_y(lMax, k_cart, E_0s[ki])
            pq_0[:,0,:] = np.exp(1j*np.sum(k_cart[ň,:]*cart_lattice,axis=-1))[:, ň] * p_y0[ň, :]
            pq_0[:,1,:] = np.exp(1j*np.sum(k_cart[ň,:]*cart_lattice,axis=-1))[:, ň] * q_y0[ň, :]
            if (return_pq_0):
                pq_0_arr[ki] = pq_0
            MP_0 = np.empty((N,2,nelem),dtype=np.complex_)
            for j in range(N): # I wonder how this can be done without this loop...
                MP_0[j] = np.tensordot(TMatrices[xij, yij],pq_0[j],axes=([-2,-1],[-2,-1]))
            MP_0.shape = (N*2*nelem,)

            ab[ki] = scipy.linalg.lu_solve(lupiv, MP_0)
        ab.shape = (K, xN, yN, 2, nelem)
    if not (return_pq_0 + return_pq + return_xy):
        return ab
    returnlist = [ab]
    if (return_pq_0):
        returnlist.append(pq_0_arr)
    if (return_pq):
        warnings.warn("return_pq not implemented, ignoring")
        # returnlist.append(pq_arr)
    if (return_xy):
        returnlist.append(x)
        returnlist.append(y)
    return tuple(returnlist) 
