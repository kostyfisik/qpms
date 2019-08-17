# [Xu] = Journal of computational physics 139, 137–165

from __future__ import print_function

def p_q(q, n, nu):
    return n + nu - 2*q

def qmax(M, n, mu, nu):
    return floor(min(n, nu, (n+nu-abs(M+mu))/2))

def Qmax(M, n, mu, nu): # [Xu](60)
    return floor(min(n, nu, (n+nu+1-abs(M+mu))/2))

def gaunta_p(M, n, mu, nu, p): # [Xu](5)
    #print (M,n,mu,nu,p, file=sys.stderr)
    return (-1)**(M+mu) * (2*p +1) * sqrt(
        factorial(n+M) * factorial(nu+mu) * factorial(p-M-mu)
      / factorial(n-M) / factorial(nu-mu) / factorial(p+M+mu)) * (
        wigner_3j(n, nu, p, 0, 0, 0) * wigner_3j(n, nu, p, M, mu, -M-mu))

def bCXcoeff(M, n, mu, nu, p): # [Xu](61)
    #print(M,n,mu,nu,p,file=sys.stderr)
    return (-1)**(M+mu) * (2*p + 3) * sqrt(
        factorial(n+M) * factorial(nu+mu) * factorial(p+1-M-mu)
      / factorial(n-M) / factorial(nu-mu) / factorial(p+1+M+mu)) * (
        wigner_3j(n, nu, p, 0, 0, 0) * wigner_3j(n, nu, p+1, M, mu, -M-mu))

def ACXcoeff(m, n, mu, nu, q): # [Xu](58)
    p = p_q(q,n,nu)
    return ((-1)**m * (2*nu + 1) * factorial(n+m) * factorial(nu-mu) / (
        2 * n * (n+1) * factorial(n-m) * factorial(nu+mu)) * I**p * 
        (n*(n+1) + nu*(nu+1) - p*(p+1)) * gaunta_p(-m,n,mu,nu,p))

def BCXcoeff(m, n, mu, nu, q): # [Xu](59)
    p = p_q(q,n,nu)
    return ((-1)**(m+1) * (2*nu + 1) * factorial(n+m) * factorial(nu-mu) / (
        2 * n * (n+1) * factorial(n-m) * factorial(nu+mu)) * I**(p+1) * 
        sqrt(((p+1)**2-(n-nu)**2) * ((n+nu+1)**2-(p+1)**2)) 
        * bCXcoeff(-m,n,mu,nu,p))

def printACXcoeffs(lMax, file=sys.stdout):
    for n in IntegerRange(lMax+1):
        for nu in IntegerRange(lMax+1):
            for m in IntegerRange(-n, n+1):
                for mu in IntegerRange(-nu, nu+1):
                    for q in IntegerRange(qmax(-m,n,mu,nu)):
                        #print(m, n, mu, nu, q, p_q(q,n,nu), file=sys.stderr)
                        coeff= ACXcoeff(m, n, mu, nu, q);
                        print(N(coeff, prec=53), 
                            ", // %d, %d, %d, %d, %d," % (m,n,mu,nu,q),
                            coeff,
                            file=file)
    return

def printBCXcoeffs(lMax, file=sys.stdout):
    for n in IntegerRange(lMax+1):
        for nu in IntegerRange(lMax+1):
            for m in IntegerRange(-n, n+1):
                for mu in IntegerRange(-nu, nu+1):
                    for q in IntegerRange(1, Qmax(-m,n,mu,nu) +1 ):
                        #print(m, n, mu, nu, q, p_q(q,n,nu), file=sys.stderr)
                        coeff= BCXcoeff(m, n, mu, nu, q);
                        print(N(coeff, prec=53), 
                            ", // %d, %d, %d, %d, %d," % (m,n,mu,nu,q),
                            coeff,
                            file=file)
    return

sphericalBessels = (None,
        spherical_bessel_J,
        spherical_bessel_Y,
        spherical_hankel1,
        spherical_hankel2
        )

# N.B. sage's gen_legendre_P _does_ include (-1)**m Condon-Shortley phase
# whereas formulae in [Xu] do not.
def trcoeff_ACX(m, n, mu, nu, besseltype, kd, th, fi, csphase=1): # [Xu](58)
    res = 0
    for q in range(qmax(-m,n,mu,nu)+1):
        p = p_q(q,n,nu)
        res += ACXcoeff(m,n,mu,nu,q) * sphericalBessels[besseltype](p,kd) * gen_legendre_P(p, mu-m, cos(th)) * (-csphase)**(mu-m) # compensate for csphase
    res *= exp(I*(mu-m)*fi)
    return res

def trcoeff_BCX(m, n, mu, nu, besseltype, kd, th, fi, csphase=1): # [Xu](59) 
    res = 0
    for q in IntegerRange(1,Qmax(-m,n,mu,nu)+1):
        p = p_q(q,n,nu)
        res += BCXcoeff(m,n,mu,nu,q) * sphericalBessels[besseltype](p+1,kd) * gen_legendre_P(p+1, mu-m, cos(th)) * (-csphase)**(mu-m)
    res *= exp(I*(mu-m)*fi)
    return res


def legpi_xu(n, m, fi): # momentálně neošetřeny okraje (cos(fi) == +- 1)
    return m/sin(fi) * gen_legendre_P(n, m, cos(fi))

def legtau_xu(n, m, fi):
    locx = var('locx')
    return -sin(fi)*derivative(gen_legendre_P(n,m,locx), locx).substitute(locx = cos(fi))

def vswf_M_xu(besseltype, n, m, kr, th, fi):
    postpart = sphericalBessels[besseltype](n, kr) * exp(I * m * fi)
    tc = I*legpi_xu(n,m,th) * postpart
    fc = -legtau_xu(n,m,th) * postpart
    return (0, tc, fc)

def vswf_N_xu(besseltype, n, m, kr, th, fi):
    eimf = exp(I * m * fi)
    rc = n * (n+1) * gen_legendre_P(n, m, cos(th)) * sphericalBessels[besseltype](n, kr)/kr * eimf
    krv = var('krv')
    radpart = derivative(krv * sphericalBessels[besseltype](n, krv), krv).substitute(krv=kr)/kr
    tc = legtau_xu(n,m,th) * radpart * eimf
    fc = I*legpi_xu(n,m,th)  * radpart * eimf
    return (rc,tc,fc)

def cart2sph(v):
    (x, y, z) = v
    r = sqrt(x**2 + y**2 + z**2)
    th = arccos(z/r) if r else 0
    fi = arctan2(y,x)
    return (r, th, fi)

def sph2cart(s):
    (r, th, fi) = s
    sinth = sin(th)
    x = r * sinth * cos(fi)
    y = r * sinth * sin(fi)
    z = r * cos(th)
    return (x,y,z)

def sphvec2cart(loccart, sph):
    r, th, fi = sph
    sinth = sin(th)
    costh = cos(th)
    sinfi = sin(fi)
    cosfi = cos(fi)

    rx = sinth * cosfi
    ry = sinth * sinfi
    rz = costh
    tx = costh * cosfi
    ty = costh * sinfi
    tz = -sinth
    fx = -sinfi
    fy = cosfi
    fz = 0
    
    rc, tc, fc = loccart
    x = rx * rc + tx * tc + fx * fc
    y = ry * rc + ty * tc + fy * fc
    z = rz * rc + tz * tc + fz * fc
    return (x, y, z)

def cart2sphvec(cart, sph):
    _, th, fi = sph
    x, y, z = cart

    rx = sinth * cosfi
    ry = sinth * sinfi
    rz = costh
    tx = costh * cosfi
    ty = costh * sinfi
    tz = -sinth
    fx = -sinfi
    fy = cosfi
    fz = 0

    rc = rx * x + ry * y + rz * z
    tc = tx * x + ty * y + tz * z
    fc = fx * x + fy * y + fz * z
    return (rc, tc, fc)

def test_M_translation_xu(lMax, origl, origm, origcartat, cartshift):
    ox, oy, oz = origcartat
    sx, sy, sz = cartshift
    newcartat = (ox - sx, oy - sy, oz - sz)
    w1s = cart2sph(origcartat)
    w2s = cart2sph(newcartat)
    pass # TODO

    

