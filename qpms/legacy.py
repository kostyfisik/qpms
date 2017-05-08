import math
import numpy as np
nx = None

_s3 = math.sqrt(3)

from qpms_c import get_mn_y, trans_calculator
from .qpms_p import cart2sph




# Functions copied from scattering_xu, additionaly normalized
from py_gmm.gmm import vec_trans as vc

#@ujit
def q_max(m,n,μ,ν):
    return min(n,ν,(n+ν-abs(m+μ))/2)

q_max_v = np.vectorize(q_max)

# returns array with indices corresponding to q
# argument q does nothing for now
#@ujit
def a_q(m,n,μ,ν,q = None):
    qm=q_max(m,n,μ,ν)
    res, err= vc.gaunt_xu(m,n,μ,ν,qm)
    if(err):
        print("m,n,μ,ν,qm = ",m,n,μ,ν,qm)
        raise ValueError('Something bad in the fortran subroutine gaunt_xu happened')
    return res

a_q_v = np.vectorize(a_q)


# All arguments are single numbers (for now)
# ZDE VYCHÁZEJÍ DIVNÁ ZNAMÉNKA
#@ujit
def Ã(m,n,μ,ν,kdlj,θlj,φlj,r_ge_d,J):
    """
    The Ã translation coefficient for spherical vector waves.
    
    Parameters
    ----------
    m, n: int
        The indices (degree and order) of the destination basis.
    μ, ν: int
        The indices of the source basis wave.
    kdlj, θlj, φlj: float
        The spherical coordinates of the relative position of
        the new center vs. the old one (R_new - R_old);
        the distance has to be already multiplied by the wavenumber!
    r_ge_d: TODO
    J: 1, 2, 3 or 4
        Type of the wave in the old center.

    Returns
    -------
    TODO

    Bugs
    ----
    gevero's gaunt coefficient implementation fails for large m, n (the unsafe territory
    is somewhere around -72, 80)

    """
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
#@ujit
def B̃(m,n,μ,ν,kdlj,θlj,φlj,r_ge_d,J):
    """
    The B̃ translation coefficient for spherical vector waves.
    
    Parameters
    ----------
    m, n: int
        The indices (degree and order) of the destination basis.
    μ, ν: int
        The indices of the source basis wave.
    kdlj, θlj, φlj: float
        The spherical coordinates of the relative position of
        the new center vs. the old one (R_new - R_old);
        the distance has to be already multiplied by the wavenumber!
    r_ge_d: TODO
    J: 1, 2, 3 or 4
        Type of the wave in the old center.

    Returns:
    --------
    TODO
    """
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




#@ujit
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
    
#@ujit
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




####################
# Array simulations
####################



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
#@ujit
def scatter_plane_wave_rectarray(omega, epsilon_b, xN, yN, xd, yd, TMatrices, k_dirs, E_0s, 
        return_pq_0 = False, return_pq= False, return_xy = False, watch_time = False):
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
    watch_time : bool
        Inform about the progress on stderr
        
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
    if (watch_time):
        timec = time.time()
        print('%.4f: running scatter_plane_wave_rectarray' % timec, file = sys.stderr)
        sys.stderr.flush()
    nelem = TMatrices.shape[-1]
    if ((nelem != TMatrices.shape[-3]) or (2 != TMatrices.shape[-2]) or (2 != TMatrices.shape[-4])):
        raise ValueError('The T-matrices must be of shape (N, 2, nelem, 2, nelem) but are of shape %s' % (str(TMatrices.shape),))
    lMax = nelem2lMax(nelem)
    if not lMax:
        raise ValueError('The "nelem" dimension of T-matrix has invalid value (%d).' % nelem)
    if (watch_time):
        print('xN = %d, yN = %d, lMax = %d' % (xN, yN, lMax), file = sys.stderr)
        sys.stderr.flush()
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
    if (watch_time):
        timec = time.time()
        print('%.4f: calculating the %d translation matrix elements' % (timec, 8*nelem*nelem*xN*yN), file = sys.stderr)
        sys.stderr.flush()
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
    if (watch_time):
        timecold = timec
        timec = time.time()
        print('%4f: translation matrix elements calculated (elapsed %.2f s), filling the matrix' 
                % (timec, timec-timecold), file = sys.stderr)
        sys.stderr.flush()
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
    if (watch_time):
        timecold = timec
        timec = time.time()
        print('%4f: translation matrix filled (elapsed %.2f s), building the interaction matrix' 
                % (timec, timec-timecold), file=sys.stderr)
        sys.stderr.flush()

    # Now we solve a linear problem (1 - M T) A = M P_0 where M is the T-matrix :-)
    MT = np.empty((N,2,nelem,N,2,nelem),dtype=np.complex_)
    
    TMatrices = np.broadcast_to(TMatrices, (xN, yN, 2, nelem, 2, nelem))
    for j in range(N): # I wonder how this can be done without this loop...
        xij, yij = xyind[j]
        MT[j] = np.tensordot(TMatrices[xij, yij],transmat[j],axes=([-2,-1],[0,1]))
    MT.shape = (N*2*nelem, N*2*nelem)
    leftmatrix = np.identity(N*2*nelem) - MT
    MT = None
    if (watch_time):
        timecold = timec
        timec = time.time()
        print('%.4f: interaction matrix complete (elapsed %.2f s)' % (timec, timec-timecold),
                file=sys.stderr)
        sys.stderr.flush()

    if ((1 == k_dirs.ndim) and (1 == E_0s.ndim)):
        k_cart = k_dirs * k_out # wave vector of the incident plane wave
        pq_0 = np.zeros((N,2,nelem), dtype=np.complex_)
        p_y0, q_y0 = plane_pq_y(lMax, k_cart, E_0s)
        pq_0[:,0,:] = np.exp(1j*np.sum(k_cart[ň,:]*cart_lattice,axis=-1))[:, ň] * p_y0[ň, :]
        pq_0[:,1,:] = np.exp(1j*np.sum(k_cart[ň,:]*cart_lattice,axis=-1))[:, ň] * q_y0[ň, :]
        if (return_pq_0):
            pq_0_arr = pq_0
        MP_0 = np.empty((N,2,nelem),dtype=np.complex_)
        #if (watch_time):
        #    print('%4f: building the interaction matrix' % time.time(), file=sys.stderr)

        for j in range(N): # I wonder how this can be done without this loop...
            MP_0[j] = np.tensordot(TMatrices[xij, yij],pq_0[j],axes=([-2,-1],[-2,-1]))
        MP_0.shape = (N*2*nelem,)
        
        if (watch_time):
            timecold = time.time()
            print('%4f: solving the scattering problem for single incoming wave' % timecold,
                    file = sys.stderr)
            sys.stderr.flush()
        ab = np.linalg.solve(leftmatrix, MP_0)
        if watch_time:
            timec = time.time()
            print('%4f: solved (elapsed %.2f s)' % (timec, timec-timecold), file=sys.stderr)
            sys.stderr.flush()

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
        if watch_time:
            timecold = time.time()
            print('%.4f: factorizing the interaction matrix' % timecold, file=sys.stderr)
            sys.stderr.flush()
        lupiv = scipy.linalg.lu_factor(leftmatrix, overwrite_a=True)
        leftmatrix = None
        if watch_time:
            timec = time.time()
            print('%.4f: factorization complete (elapsed %.2f s)' % (timec, timec-timecold),
                    file = sys.stderr)
            print('%.4f: solving the scattering problem for %d incoming waves' % (timec, K),
                    file=sys.stderr)
            sys.stderr.flush()
            timecold = timec
        
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
        if watch_time:
            timec = time.time()
            print('%.4f: done (elapsed %.2f s)' % (timec, timec-timecold),file = sys.stderr)
            sys.stderr.flush()
    if not (return_pq_0 + return_pq + return_xy):
        return ab
    returnlist = [ab]
    if (return_pq_0):
        pq_0_arr.shape = ab.shape
        returnlist.append(pq_0_arr)
    if (return_pq):
        warnings.warn("return_pq not implemented, ignoring")
        # returnlist.append(pq_arr)
    if (return_xy):
        returnlist.append(x)
        returnlist.append(y)
    return tuple(returnlist) 


import warnings
#@ujit
def scatter_constmultipole_rectarray(omega, epsilon_b, xN, yN, xd, yd, TMatrices, pq_0_c = 1, 
        return_pq= False, return_xy = False, watch_time = False):
    """
    Solves the plane wave linear scattering problem for a rectangular array of particles
    for one frequency and constant exciting spherical waves throughout the array.
    
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
    pq_0_c : (nelem)-shaped complex array or compatible
        The initial excitation coefficients for the ("complex") multipole waves, in Taylor's convention.
    return_pq : bool NOT IMPLEMENTED
        Return also the multipole decomposition coefficients of the field incoming to each
        particle (inc. the field scattered from other particles.
    return_xy : bool
        Return also the cartesian x, y positions of the particles. 
    watch_time : bool
        Inform about the progress on stderr
        
    Returns
    -------
    ab : (nelem, xN, yN, 2, nelem)-shaped complex array
        The a (electric wave), b (magnetic wave) coefficients of the outgoing field for each particle.
        If none of return_pq or return_xy is set, the array is not enclosed in a tuple.
    pq : (nelem, xN, yN, 2, nelem)-shaped complex array NOT IMPLEMENTED
        The p (electric wave), q (magnetic wave) coefficients of the total exciting field
        for each particle (including the field scattered from other particles)
    x, y : (xN, yN)-shaped real array
        The x,y positions of the nanoparticles.
    """
    if (watch_time):
        timec = time.time()
        print('%.4f: running scatter_plane_wave_rectarray' % timec, file = sys.stderr)
        sys.stderr.flush()
    nelem = TMatrices.shape[-1]
    if ((nelem != TMatrices.shape[-3]) or (2 != TMatrices.shape[-2]) or (2 != TMatrices.shape[-4])):
        raise ValueError('The T-matrices must be of shape (N, 2, nelem, 2, nelem) but are of shape %s' % (str(TMatrices.shape),))
    lMax = nelem2lMax(nelem)
    if not lMax:
        raise ValueError('The "nelem" dimension of T-matrix has invalid value (%d).' % nelem)
    if (watch_time):
        print('xN = %d, yN = %d, lMax = %d' % (xN, yN, lMax), file = sys.stderr)
        sys.stderr.flush()
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
    if (watch_time):
        timec = time.time()
        print('%.4f: calculating the %d translation matrix elements' % (timec, 8*nelem*nelem*xN*yN), file = sys.stderr)
        sys.stderr.flush()
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
    if (watch_time):
        timecold = timec
        timec = time.time()
        print('%4f: translation matrix elements calculated (elapsed %.2f s), filling the matrix' 
                % (timec, timec-timecold), file = sys.stderr)
        sys.stderr.flush()
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
    if (watch_time):
        timecold = timec
        timec = time.time()
        print('%4f: translation matrix filled (elapsed %.2f s), building the interaction matrix' 
                % (timec, timec-timecold), file=sys.stderr)
        sys.stderr.flush()

    # Now we solve a linear problem (1 - M T) A = M P_0 where M is the T-matrix :-)
    MT = np.empty((N,2,nelem,N,2,nelem),dtype=np.complex_)
    
    TMatrices = np.broadcast_to(TMatrices, (xN, yN, 2, nelem, 2, nelem))
    for j in range(N): # I wonder how this can be done without this loop...
        xij, yij = xyind[j]
        MT[j] = np.tensordot(TMatrices[xij, yij],transmat[j],axes=([-2,-1],[0,1]))
    MT.shape = (N*2*nelem, N*2*nelem)
    leftmatrix = np.identity(N*2*nelem) - MT
    MT = None
    if (watch_time):
        timecold = timec
        timec = time.time()
        print('%.4f: interaction matrix complete (elapsed %.2f s)' % (timec, timec-timecold),
                file=sys.stderr)
        sys.stderr.flush()

    # А ну, чики-брики и в дамки!
    if watch_time:
        timecold = time.time()
        print('%.4f: factorizing the interaction matrix' % timecold, file=sys.stderr)
        sys.stderr.flush()
    lupiv = scipy.linalg.lu_factor(leftmatrix, overwrite_a=True)
    leftmatrix = None
    if watch_time:
        timec = time.time()
        print('%.4f: factorization complete (elapsed %.2f s)' % (timec, timec-timecold),
                file = sys.stderr)
        print('%.4f: solving the scattering problem for %d incoming multipoles' % (timec, nelem*2),
                file=sys.stderr)
        sys.stderr.flush()
        timecold = timec
    
    if(pq_0_c == 1):
        pq_0_c = np.full((2,nelem),1)
    ab = np.empty((2,nelem,N*2*nelem), dtype=complex)
    for N_or_M in range(2):
      for yy in range(nelem):
        pq_0 = np.zeros((2,nelem), dtype=np.complex_)
        pq_0[N_or_M,yy] = pq_0_c[N_or_M,yy]
        pq_0 = np.broadcast_to(pq_0, (N, 2, nelem))
        MP_0 = np.empty((N,2,nelem),dtype=np.complex_)
        for j in range(N): # I wonder how this can be done without this loop...
            xij, yij = xyind[j]
            MP_0[j] = np.tensordot(TMatrices[xij, yij],pq_0[j],axes=([-2,-1],[-2,-1]))
        MP_0.shape = (N*2*nelem,)

        ab[N_or_M, yy] = scipy.linalg.lu_solve(lupiv, MP_0)
    ab.shape = (2,nelem, xN, yN, 2, nelem)
    if watch_time:
        timec = time.time()
        print('%.4f: done (elapsed %.2f s)' % (timec, timec-timecold),file = sys.stderr)
        sys.stderr.flush()
    if not (return_pq + return_xy):
        return ab
    returnlist = [ab]
    if (return_pq):
        warnings.warn("return_pq not implemented, ignoring")
        # returnlist.append(pq_arr)
    if (return_xy):
        returnlist.append(x)
        returnlist.append(y)
    return tuple(returnlist) 













# --------------    hexagonal lattice translation coefficients   ---------------------------
# Implementation using (C) one-by-one AB-coefficient calculation ufunc
def hexlattice_precalc_AB_save2(file, lMax, k_hexside, maxlayer, circular=True, savepointinfo = False, J_scat=3):
    params = {
        'lMax' : lMax,
        'k_hexside' : k_hexside,
        'maxlayer' : maxlayer,
        'circular' : circular,
        'savepointinfo' : savepointinfo,
        'J_scat' : J_scat
    }
    tpdict = generate_trianglepoints(maxlayer, v3d=True, circular=circular, sixthindices=True, mirrorindices=True)
    tphcdict = generate_trianglepoints_hexcomplement(maxlayer, v3d=True, circular=circular, thirdindices=True, mirrorindices=True)
    my, ny = get_mn_y(lMax)
    nelem = len(my)
    a_self_nm = np.empty((tpdict['nmi'].shape[0],nelem,nelem), dtype=complex)
    b_self_nm = np.empty((tpdict['nmi'].shape[0],nelem,nelem), dtype=complex)
    a_self_m0 = np.empty((tpdict['mi'].shape[1],nelem,nelem), dtype=complex)
    b_self_m0 = np.empty((tpdict['mi'].shape[1],nelem,nelem), dtype=complex)
    a_d2u_nm = np.empty((tphcdict['nmi'].shape[0],nelem,nelem), dtype=complex)
    b_d2u_nm = np.empty((tphcdict['nmi'].shape[0],nelem,nelem), dtype=complex)
    a_d2u_m0 = np.empty((tphcdict['mi'].shape[1],nelem,nelem), dtype=complex)
    b_d2u_m0 = np.empty((tphcdict['mi'].shape[1],nelem,nelem), dtype=complex)
    
    k_0 = k_hexside*_s3 # not really a wave vector here because of the normalisation!
    tc = trans_calculator(lMax)
    
    y = np.arange(nelem)

    points = tpdict['points'][tpdict['nmi']]
    d_i2j = cart2sph(points)
    a_self_nm, b_self_nm = tc.get_AB(my[nx,:,nx],ny[nx,:,nx],my[nx,nx,:],ny[nx,nx,:],k_0*d_i2j[:,nx,nx,0],d_i2j[:,nx,nx,1],d_i2j[:,nx,nx,2],False,J_scat)

    points = tpdict['points'][tpdict['mi'][0]]
    d_i2j = cart2sph(points)
    a_self_m0, b_self_m0 = tc.get_AB(my[nx,:,nx],ny[nx,:,nx],my[nx,nx,:],ny[nx,nx,:],k_0*d_i2j[:,nx,nx,0],d_i2j[:,nx,nx,1],d_i2j[:,nx,nx,2],False,J_scat)
   
    points = tphcdict['points'][tphcdict['nmi']]
    d_i2j = cart2sph(points)
    a_d2u_nm, b_d2u_nm = tc.get_AB(my[nx,:,nx],ny[nx,:,nx],my[nx,nx,:],ny[nx,nx,:],k_0*d_i2j[:,nx,nx,0],d_i2j[:,nx,nx,1],d_i2j[:,nx,nx,2],False,J_scat)
  
    points = tphcdict['points'][tphcdict['mi'][0]]
    d_i2j = cart2sph(points)
    a_d2u_m0, b_d2u_m0 = tc.get_AB(my[nx,:,nx],ny[nx,:,nx],my[nx,nx,:],ny[nx,nx,:],k_0*d_i2j[:,nx,nx,0],d_i2j[:,nx,nx,1],d_i2j[:,nx,nx,2],False,J_scat)
  
    tosave = {
        'a_self_nm' : a_self_nm,
        'a_self_m0' : a_self_m0,
        'b_self_nm' : b_self_nm,
        'b_self_m0' : b_self_m0,
        'a_d2u_nm' : a_d2u_nm,
        'a_d2u_m0' : a_d2u_m0,
        'b_d2u_nm' : b_d2u_nm,
        'b_d2u_m0' : b_d2u_m0,
        'precalc_params' : params
    }
    if savepointinfo:
        tosave['tp_points'] = tpdict['points'],
        tosave['tp_si'] = tpdict['si'],
        tosave['tp_mi'] = tpdict['mi'],
        tosave['tp_nmi'] = tpdict['nmi']
        tosave['tphc_points'] = tphcdict['points'],
        tosave['tphc_ti'] = tphcdict['ti'],
        tosave['tphc_mi'] = tphcdict['mi'],
        tosave['tphc_nmi'] = tphcdict['nmi']
    np.savez(file, **tosave)



# The oldest implementation, using the super-inefficient pure python translation coefficients
def hexlattice_precalc_AB_save_purepy(file, lMax, k_hexside, maxlayer, circular=True, savepointinfo = False, J_scat=3):
    params = {
        'lMax' : lMax,
        'k_hexside' : k_hexside,
        'maxlayer' : maxlayer,
        'circular' : circular,
        'savepointinfo' : savepointinfo,
        'J_scat' : J_scat
    }
    tpdict = generate_trianglepoints(maxlayer, v3d=True, circular=circular, sixthindices=True, mirrorindices=True)
    tphcdict = generate_trianglepoints_hexcomplement(maxlayer, v3d=True, circular=circular, thirdindices=True, mirrorindices=True)
    my, ny = get_mn_y(lMax)
    nelem = len(my)
    a_self_nm = np.empty((tpdict['nmi'].shape[0],nelem,nelem), dtype=complex)
    b_self_nm = np.empty((tpdict['nmi'].shape[0],nelem,nelem), dtype=complex)
    a_self_m0 = np.empty((tpdict['mi'].shape[1],nelem,nelem), dtype=complex)
    b_self_m0 = np.empty((tpdict['mi'].shape[1],nelem,nelem), dtype=complex)
    a_d2u_nm = np.empty((tphcdict['nmi'].shape[0],nelem,nelem), dtype=complex)
    b_d2u_nm = np.empty((tphcdict['nmi'].shape[0],nelem,nelem), dtype=complex)
    a_d2u_m0 = np.empty((tphcdict['mi'].shape[1],nelem,nelem), dtype=complex)
    b_d2u_m0 = np.empty((tphcdict['mi'].shape[1],nelem,nelem), dtype=complex)
    
    k_0 = k_hexside*_s3 # not really a wave vector here because of the normalisation!
    
    points = tpdict['points'][tpdict['nmi']]
    for j in range(points.shape[0]):
        d_i2j = cart2sph(points[j])
        for yi in range(nelem):
            for yj in range(nelem):
                a_self_nm[j, yj, yi] = Ã(my[yj],ny[yj],my[yi],ny[yi],kdlj=d_i2j[0]*k_0,θlj=d_i2j[1],φlj=d_i2j[2],r_ge_d=False,J=J_scat)
                b_self_nm[j, yj, yi] = B̃(my[yj],ny[yj],my[yi],ny[yi],kdlj=d_i2j[0]*k_0,θlj=d_i2j[1],φlj=d_i2j[2],r_ge_d=False,J=J_scat)
    points = tpdict['points'][tpdict['mi'][0]]
    for j in range(points.shape[0]):
        d_i2j = cart2sph(points[j])
        for yi in range(nelem):
            for yj in range(nelem):
                a_self_m0[j, yj, yi] = Ã(my[yj],ny[yj],my[yi],ny[yi],kdlj=d_i2j[0]*k_0,θlj=d_i2j[1],φlj=d_i2j[2],r_ge_d=False,J=J_scat)
                b_self_m0[j, yj, yi] = B̃(my[yj],ny[yj],my[yi],ny[yi],kdlj=d_i2j[0]*k_0,θlj=d_i2j[1],φlj=d_i2j[2],r_ge_d=False,J=J_scat)
    points = tphcdict['points'][tphcdict['nmi']]
    for j in range(points.shape[0]):
        d_i2j = cart2sph(points[j])
        for yi in range(nelem):
            for yj in range(nelem):
                a_d2u_nm[j, yj, yi] = Ã(my[yj],ny[yj],my[yi],ny[yi],kdlj=d_i2j[0]*k_0,θlj=d_i2j[1],φlj=d_i2j[2],r_ge_d=False,J=J_scat)
                b_d2u_nm[j, yj, yi] = B̃(my[yj],ny[yj],my[yi],ny[yi],kdlj=d_i2j[0]*k_0,θlj=d_i2j[1],φlj=d_i2j[2],r_ge_d=False,J=J_scat)
    points = tphcdict['points'][tphcdict['mi'][0]]
    for j in range(points.shape[0]):
        d_i2j = cart2sph(points[j])
        for yi in range(nelem):
            for yj in range(nelem):
                a_d2u_m0[j, yj, yi] = Ã(my[yj],ny[yj],my[yi],ny[yi],kdlj=d_i2j[0]*k_0,θlj=d_i2j[1],φlj=d_i2j[2],r_ge_d=False,J=J_scat)
                b_d2u_m0[j, yj, yi] = B̃(my[yj],ny[yj],my[yi],ny[yi],kdlj=d_i2j[0]*k_0,θlj=d_i2j[1],φlj=d_i2j[2],r_ge_d=False,J=J_scat)
    tosave = {
        'a_self_nm' : a_self_nm,
        'a_self_m0' : a_self_m0,
        'b_self_nm' : b_self_nm,
        'b_self_m0' : b_self_m0,
        'a_d2u_nm' : a_d2u_nm,
        'a_d2u_m0' : a_d2u_m0,
        'b_d2u_nm' : b_d2u_nm,
        'b_d2u_m0' : b_d2u_m0,
        'precalc_params' : params
    }
    if savepointinfo:
        tosave['tp_points'] = tpdict['points'],
        tosave['tp_si'] = tpdict['si'],
        tosave['tp_mi'] = tpdict['mi'],
        tosave['tp_nmi'] = tpdict['nmi']
        tosave['tphc_points'] = tphcdict['points'],
        tosave['tphc_ti'] = tphcdict['ti'],
        tosave['tphc_mi'] = tphcdict['mi'],
        tosave['tphc_nmi'] = tphcdict['nmi']
    np.savez(file, **tosave)

def hexlattice_precalc_AB_loadunwrap(file, tpdict = None, tphcdict = None, return_points = False):
    npz = np.load(file)
    precalc_params = npz['precalc_params'][()]
    my, ny = get_mn_y(precalc_params['lMax'])
    nelem = len(my)
    # this I should have made more universal...
    if precalc_params['savepointinfo']:
        if not tpdict:
            tpdict = {
                'points' : npz['tp_points'],
                'si' : npz['tp_si'],
                'mi' : npz['tp_mi'],
                'nmi' : npz['tp_nmi'],
            }
        if not tphcdict:
            tphcdict = {
                'points' : npz['tphc_points'],
                'ti' : npz['tphc_ti'],
                'mi' : npz['tphc_mi'],
                'nmi' : npz['tphc_nmi']            
        }
    else:
        if not tpdict:
            tpdict = generate_trianglepoints(maxlayer = precalc_params['maxlayer'], v3d=True, 
                                             circular=precalc_params['circular'], sixthindices=True, mirrorindices=True)
        if not tphcdict:
            tphcdict = generate_trianglepoints_hexcomplement(maxlayer=precalc_params['maxlayer'], v3d=True, 
                                                             circular=precalc_params['circular'], thirdindices=True, mirrorindices=True)
        
    # For some obscure reason, I keep getting trailing single-dimension in the beginning for these arrays
    for a in (tpdict['points'], tphcdict['points'], tpdict['si'], tpdict['mi'],
             tphcdict['ti'], tphcdict['mi']):
        if len(a.shape) > 2:
            a.shape = a.shape[1::] 
        
    self_tr = tpdict['points'] 
    d2u_tr = tphcdict['points']
    if len(self_tr.shape)>2:
        self_tr = np.reshape(self_tr, self_tr.shape[1::])
    if len(d2u_tr.shape)>2:
        d2u_tr = np.reshape(d2u_tr, d2u_tr.shape[1::])
    u2d_tr = -d2u_tr
    a_self = np.empty((self_tr.shape[0],nelem,nelem), dtype=complex)
    b_self = np.empty((self_tr.shape[0],nelem,nelem), dtype=complex)
    a_d2u  = np.empty(( d2u_tr.shape[0],nelem,nelem), dtype=complex)
    b_d2u  = np.empty(( d2u_tr.shape[0],nelem,nelem), dtype=complex)
    a_self[tpdict['nmi']]=npz['a_self_nm']
    a_self[tpdict['mi'][0]]=npz['a_self_m0']
    b_self[tpdict['nmi']]=npz['b_self_nm']
    b_self[tpdict['mi'][0]]=npz['b_self_m0']
    mirrorangles = cart2sph(self_tr[tpdict['mi'][1]])[:,2] - cart2sph(self_tr[tpdict['mi'][0]])[:,2]
    a_self[tpdict['mi'][1],:,:] = a_self[tpdict['mi'][0],:,:] * np.exp(1j*mirrorangles[:,nx,nx]*(my[nx,nx,:]-my[nx,:,nx]))
    b_self[tpdict['mi'][1],:,:] = b_self[tpdict['mi'][0],:,:] * np.exp(1j*mirrorangles[:,nx,nx]*(my[nx,nx,:]-my[nx,:,nx]))
    for i in range(1,6):
        a_self[tpdict['si'][i],:,:] = a_self[tpdict['si'][0],:,:] * np.exp(1j*math.pi/3*i*(my[nx,:]-my[:,nx]))
        b_self[tpdict['si'][i],:,:] = b_self[tpdict['si'][0],:,:] * np.exp(1j*math.pi/3*i*(my[nx,:]-my[:,nx]))
    a_d2u[tphcdict['nmi']]=npz['a_d2u_nm']
    a_d2u[tphcdict['mi'][0]]=npz['a_d2u_m0']
    b_d2u[tphcdict['nmi']]=npz['b_d2u_nm']
    b_d2u[tphcdict['mi'][0]]=npz['b_d2u_m0']
    mirrorangles = cart2sph(self_tr[tphcdict['mi'][1]])[:,2] - cart2sph(self_tr[tphcdict['mi'][0]])[:,2]
    a_d2u[tphcdict['mi'][1],:,:] = a_d2u[tphcdict['mi'][0],:,:] * np.exp(1j*mirrorangles[:,nx,nx]*(my[nx,nx,:]-my[nx,:,nx]))
    b_d2u[tphcdict['mi'][1],:,:] = b_d2u[tphcdict['mi'][0],:,:] * np.exp(1j*mirrorangles[:,nx,nx]*(my[nx,nx,:]-my[nx,:,nx]))
    for i in (1,-1):
        a_d2u[tphcdict['ti'][i],:,:] = a_d2u[tphcdict['ti'][0],:,:] * np.exp(i*2j*math.pi/3*(my[nx,:]-my[:,nx]))
        b_d2u[tphcdict['ti'][i],:,:] = b_d2u[tphcdict['ti'][0],:,:] * np.exp(i*2j*math.pi/3*(my[nx,:]-my[:,nx]))
    a_u2d = a_d2u * (-1)**(my[nx,:]-my[:,nx])
    b_u2d = b_d2u * (-1)**(my[nx,:]-my[:,nx])
    d = {
        'a_self' : a_self,
        'b_self' : b_self,
        'a_d2u'  : a_d2u,
        'b_d2u'  : b_d2u,
        'a_u2d'  : a_u2d,
        'b_u2d'  : b_u2d,
    }
    for k in precalc_params.keys():
        d[k] = precalc_params[k]
    if return_points:
        d['d2u_tr'] = tphcdict['points']
        d['u2d_tr'] = -tphcdict['points']
        d['self_tr'] = tpdict['points']
    return d
    
def hexlattice_get_AB(lMax, k_hexside, maxlayer, circular=True, return_points = True, J_scat=3):
    params = {
        'lMax' : lMax,
        'k_hexside' : k_hexside,
        'maxlayer' : maxlayer,
        'circular' : circular,
        'savepointinfo' : return_points, # should I delete this key?
        'J_scat' : J_scat
    }
    tpdict = generate_trianglepoints(maxlayer, v3d=True, circular=circular, sixthindices=True, mirrorindices=True)
    tphcdict = generate_trianglepoints_hexcomplement(maxlayer, v3d=True, circular=circular, thirdindices=True, mirrorindices=True)
    my, ny = get_mn_y(lMax)
    nelem = len(my)
    a_self_nm = np.empty((tpdict['nmi'].shape[0],nelem,nelem), dtype=complex)
    b_self_nm = np.empty((tpdict['nmi'].shape[0],nelem,nelem), dtype=complex)
    a_self_m0 = np.empty((tpdict['mi'].shape[1],nelem,nelem), dtype=complex)
    b_self_m0 = np.empty((tpdict['mi'].shape[1],nelem,nelem), dtype=complex)
    a_d2u_nm = np.empty((tphcdict['nmi'].shape[0],nelem,nelem), dtype=complex)
    b_d2u_nm = np.empty((tphcdict['nmi'].shape[0],nelem,nelem), dtype=complex)
    a_d2u_m0 = np.empty((tphcdict['mi'].shape[1],nelem,nelem), dtype=complex)
    b_d2u_m0 = np.empty((tphcdict['mi'].shape[1],nelem,nelem), dtype=complex)
    
    k_0 = k_hexside*_s3 # not really a wave vector here because of the normalisation!
    tc = trans_calculator(lMax)
    
    y = np.arange(nelem)

    points = tpdict['points'][tpdict['nmi']]
    d_i2j = cart2sph(points)
    a_self_nm, b_self_nm = tc.get_AB_arrays(k_0*d_i2j[:,0],d_i2j[:,1],d_i2j[:,2],np.array([False]),J_scat)

    points = tpdict['points'][tpdict['mi'][0]]
    d_i2j = cart2sph(points)
    a_self_m0, b_self_m0 = tc.get_AB_arrays(k_0*d_i2j[:,0],d_i2j[:,1],d_i2j[:,2],np.array([False]),J_scat)
   
    points = tphcdict['points'][tphcdict['nmi']]
    d_i2j = cart2sph(points)
    a_d2u_nm, b_d2u_nm = tc.get_AB_arrays(k_0*d_i2j[:,0],d_i2j[:,1],d_i2j[:,2],np.array([False]),J_scat)
  
    points = tphcdict['points'][tphcdict['mi'][0]]
    d_i2j = cart2sph(points)
    a_d2u_m0, b_d2u_m0 = tc.get_AB_arrays(k_0*d_i2j[:,0],d_i2j[:,1],d_i2j[:,2],np.array([False]),J_scat)
    '''
    tosave = {
        'a_self_nm' : a_self_nm,
        'a_self_m0' : a_self_m0,
        'b_self_nm' : b_self_nm,
        'b_self_m0' : b_self_m0,
        'a_d2u_nm' : a_d2u_nm,
        'a_d2u_m0' : a_d2u_m0,
        'b_d2u_nm' : b_d2u_nm,
        'b_d2u_m0' : b_d2u_m0,
        'precalc_params' : params
    }
    if savepointinfo:
        tosave['tp_points'] = tpdict['points'],
        tosave['tp_si'] = tpdict['si'],
        tosave['tp_mi'] = tpdict['mi'],
        tosave['tp_nmi'] = tpdict['nmi']
        tosave['tphc_points'] = tphcdict['points'],
        tosave['tphc_ti'] = tphcdict['ti'],
        tosave['tphc_mi'] = tphcdict['mi'],
        tosave['tphc_nmi'] = tphcdict['nmi']
    np.savez(file, **tosave)
    '''
    self_tr = tpdict['points'] 
    d2u_tr = tphcdict['points']
    if len(self_tr.shape)>2:
        self_tr = np.reshape(self_tr, self_tr.shape[1::])
    if len(d2u_tr.shape)>2:
        d2u_tr = np.reshape(d2u_tr, d2u_tr.shape[1::])
    u2d_tr = -d2u_tr
    a_self = np.empty((self_tr.shape[0],nelem,nelem), dtype=complex)
    b_self = np.empty((self_tr.shape[0],nelem,nelem), dtype=complex)
    a_d2u  = np.empty(( d2u_tr.shape[0],nelem,nelem), dtype=complex)
    b_d2u  = np.empty(( d2u_tr.shape[0],nelem,nelem), dtype=complex)
    a_self[tpdict['nmi']]=a_self_nm
    a_self[tpdict['mi'][0]]=a_self_m0
    b_self[tpdict['nmi']]=b_self_nm
    b_self[tpdict['mi'][0]]=b_self_m0
    mirrorangles = cart2sph(self_tr[tpdict['mi'][1]])[:,2] - cart2sph(self_tr[tpdict['mi'][0]])[:,2]
    a_self[tpdict['mi'][1],:,:] = a_self[tpdict['mi'][0],:,:] * np.exp(1j*mirrorangles[:,nx,nx]*(my[nx,nx,:]-my[nx,:,nx]))
    b_self[tpdict['mi'][1],:,:] = b_self[tpdict['mi'][0],:,:] * np.exp(1j*mirrorangles[:,nx,nx]*(my[nx,nx,:]-my[nx,:,nx]))
    for i in range(1,6):
        a_self[tpdict['si'][i],:,:] = a_self[tpdict['si'][0],:,:] * np.exp(1j*math.pi/3*i*(my[nx,:]-my[:,nx]))
        b_self[tpdict['si'][i],:,:] = b_self[tpdict['si'][0],:,:] * np.exp(1j*math.pi/3*i*(my[nx,:]-my[:,nx]))
    a_d2u[tphcdict['nmi']]=a_d2u_nm
    a_d2u[tphcdict['mi'][0]]=a_d2u_m0
    b_d2u[tphcdict['nmi']]=b_d2u_nm
    b_d2u[tphcdict['mi'][0]]=b_d2u_m0
    mirrorangles = cart2sph(self_tr[tphcdict['mi'][1]])[:,2] - cart2sph(self_tr[tphcdict['mi'][0]])[:,2]
    a_d2u[tphcdict['mi'][1],:,:] = a_d2u[tphcdict['mi'][0],:,:] * np.exp(1j*mirrorangles[:,nx,nx]*(my[nx,nx,:]-my[nx,:,nx]))
    b_d2u[tphcdict['mi'][1],:,:] = b_d2u[tphcdict['mi'][0],:,:] * np.exp(1j*mirrorangles[:,nx,nx]*(my[nx,nx,:]-my[nx,:,nx]))
    for i in (1,-1):
        a_d2u[tphcdict['ti'][i],:,:] = a_d2u[tphcdict['ti'][0],:,:] * np.exp(i*2j*math.pi/3*(my[nx,:]-my[:,nx]))
        b_d2u[tphcdict['ti'][i],:,:] = b_d2u[tphcdict['ti'][0],:,:] * np.exp(i*2j*math.pi/3*(my[nx,:]-my[:,nx]))
    a_u2d = a_d2u * (-1)**(my[nx,:]-my[:,nx])
    b_u2d = b_d2u * (-1)**(my[nx,:]-my[:,nx])
    d = {
        'a_self' : a_self,
        'b_self' : b_self,
        'a_d2u'  : a_d2u,
        'b_d2u'  : b_d2u,
        'a_u2d'  : a_u2d,
        'b_u2d'  : b_u2d,
    }
    for k in params.keys():
        d[k] = params[k]
    if return_points:
        d['d2u_tr'] = tphcdict['points']
        d['u2d_tr'] = -tphcdict['points']
        d['self_tr'] = tpdict['points']
    return d
 

