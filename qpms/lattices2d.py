import numpy as np
from enum import Enum
from math import floor

nx = None

class LatticeType(Enum):
    """
    All the five Bravais lattices in 2D
    """
    OBLIQUE=1
    RECTANGULAR=2
    SQUARE=4
    RHOMBIC=5
    EQUILATERAL_TRIANGULAR=3
    RIGHT_ISOSCELES=SQUARE
    PARALLELOGRAMMIC=OBLIQUE
    CENTERED_RHOMBIC=RECTANGULAR
    RIGHT_TRIANGULAR=RECTANGULAR
    CENTERED_RECTANGULAR=RHOMBIC
    ISOSCELE_TRIANGULAR=RHOMBIC
    RIGHT_ISOSCELE_TRIANGULAR=SQUARE
    HEXAGONAL=EQUILATERAL_TRIANGULAR

def reduceBasisSingle(b1, b2):
    """
    Lagrange-Gauss reduction of a 2D basis.
    cf. https://www.math.auckland.ac.nz/~sgal018/crypto-book/ch17.pdf
    inputs and outputs are (2,)-shaped numpy arrays
    The output shall satisfy |b1| <= |b2| <= |b2 - b1| 
    TODO doc
    
    TODO perhaps have the (on-demand?) guarantee of obtuse angle between b1, b2?
    TODO possibility of returning the (in-order, no-obtuse angles) b as well?
    """
    b1 = np.array(b1)
    b2 = np.array(b2)
    if b1.shape != (2,) or b2.shape != (2,):
        raise ValueError('Shape of b1 and b2 must be (2,)')
    B1 = np.sum(b1 * b1, axis=-1, keepdims=True)
    mu = np.sum(b1 * b2, axis=-1, keepdims=True) / B1
    b2 = b2 - np.rint(mu) * b1
    B2 = np.sum(b2 * b2, axis=-1, keepdims=True)
    while(np.any(B2 < B1)):
        b2t = b1
        b1 = b2
        b2 = b2t
        B1 = B2
        mu = np.sum(b1 * b2, axis=-1, keepdims=True) / B1
        b2 = b2 - np.rint(mu) * b1
        B2 = np.sum(b2*b2, axis=-1, keepdims=True)
    return np.array((b1,b2))

def shortestBase3(b1, b2):
    ''' 
    returns the "ordered shortest triple" of base vectors (each pair from 
    the triple is a base) and there may not be obtuse angle between b1, b2
    and between b2, b3
    '''
    b1, b2 = reduceBasisSingle(b1,b2)
    if is_obtuse(b1, b2, tolerance=0):
        b3 = b2
        b2 = b2 + b1
    else:
        b3 = b2 - b1
    return (b1, b2, b3)

def shortestBase46(b1, b2, tolerance=1e-13):
    b1, b2 = reduceBasisSingle(b1,b2)
    b1s = np.sum(b1 ** 2)
    b2s = np.sum(b2 ** 2)
    b3 = b2 - b1
    b3s = np.sum(b3 ** 2)
    eps = tolerance * (b2s + b1s)
    if abs(b3s - b2s - b1s) < eps:
        return(b1, b2, -b1, -b2)
    else:
        if b3s - b2s - b1s > eps: #obtuse
            b3 = b2
            b2 = b2 + b1
        return (b1, b2, b3, -b1, -b2, -b3)
    

def is_obtuse(b1, b2, tolerance=1e-13):
    b1s = np.sum(b1 ** 2)
    b2s = np.sum(b2 ** 2)
    b3 = b2 - b1
    b3s = np.sum(b3 ** 2)
    eps = tolerance * (b2s + b1s)
    return (b3s - b2s - b1s > eps)

def classifyLatticeSingle(b1, b2, tolerance=1e-13):
    """
    Given two basis vectors, returns 2D Bravais lattice type.
    Tolerance is relative.
    TODO doc
    """
    b1, b2 = reduceBasisSingle(b1, b2)
    b1s = np.sum(b1 ** 2)
    b2s = np.sum(b2 ** 2)
    b3 = b2 - b1
    b3s = np.sum(b3 ** 2)
    eps = tolerance * (b2s + b1s)
    # Avoid obtuse angle between b1 and b2. TODO This should be yet thoroughly tested.
    # TODO use is_obtuse here?
    if b3s - b2s - b1s > eps:
        b3 = b2
        b2 = b2 + b1
        # N. B. now the assumption |b3| >= |b2| is no longer valid
        #b3 = b2 - b1
        b2s = np.sum(b2 ** 2)
        b3s = np.sum(b3 ** 2)
    if abs(b2s - b1s) < eps or abs(b2s - b3s) < eps: # isoscele
        if abs(b3s - b1s) < eps:
            return LatticeType.EQUILATERAL_TRIANGULAR
        elif abs(b3s - 2 * b1s) < eps:
            return LatticeType.SQUARE
        else:
            return LatticeType.RHOMBIC
    elif abs(b3s - b2s - b1s) < eps:
        return LatticeType.RECTANGULAR
    else:
        return LatticeType.OBLIQUE

def range2D(maxN, mini=1, minj=0, minN = 0):
    """
    "Triangle indices"
    Generates pairs of non-negative integer indices (i, j) such that
    minN ≤ i + j ≤ maxN, i ≥ mini, j ≥ minj.
    TODO doc and possibly different orderings
    """
    for maxn in range(min(mini, minj, minN), floor(maxN+1)): # i + j == maxn
        for i in range(mini, maxn + 1):
            yield (i, maxn - i)

            
def generateLattice(b1, b2, maxlayer=5, include_origin=False, order='leaves'):
    bvs = shortestBase46(b1, b2)
    cc = len(bvs) # "corner count"
    
    if order == 'leaves':
        indices = np.array(list(range2D(maxlayer)))
        ia = indices[:,0]
        ib = indices[:,1]
        cc = len(bvs) # 4 for square/rec,
        leaves = list()
        if include_origin: leaves.append(np.array([[0,0]]))
        for c in range(cc):
            ba = bvs[c]
            bb = bvs[(c+1)%cc]
            leaves.append(ia[:,nx]*ba + ib[:,nx]*bb)
        return np.concatenate(leaves)
    else: 
        raise ValueError('Lattice point order not implemented: ', order)
        
def generateLatticeDisk(b1, b2, r, include_origin=False, order='leaves'):
    b1, b2 = reduceBasisSingle(b1,b2)
    blen = np.linalg.norm(b1, ord=2)
    maxlayer = 2*r/blen # FIXME kanon na vrabce? Nestačí odmocnina ze 2?
    points = generateLattice(b1,b2, maxlayer=maxlayer, include_origin=include_origin, order=order)
    mask = (np.linalg.norm(points, axis=-1, ord=2) <= r)
    return points[mask]

def cellCornersWS(b1, b2,):
    """
    Given basis vectors, returns the corners of the Wigner-Seitz unit cell
    (w1, w2, -w1, w2) for rectangular and square lattice or
    (w1, w2, w3, -w1, -w2, -w3) otherwise
    """
    def solveWS(v1, v2):
        v1x = v1[0]
        v1y = v1[1]
        v2x = v2[0]
        v2y = v2[1]
        lsm = ((-v1y, v2y), (v1x, -v2x))
        rs = ((v1x-v2x)/2, (v1y - v2y)/2)
        t = np.linalg.solve(lsm, rs)
        return np.array(v1)/2 + t[0]*np.array((v1y, -v1x))
    b1, b2 = reduceBasisSingle(b1, b2)
    latticeType = classifyLatticeSingle(b1, b2)
    if latticeType is LatticeType.RECTANGULAR or latticeType is LatticeType.SQUARE:
        return np.array( (
            (+b1+b2),
            (+b2-b1),
            (-b1-b2),
            (-b2+b1),
        )) / 2
    else:
        bvs = shortestBase46(b1,b2,tolerance=0)
        return np.array([solveWS(bvs[i], bvs[(i+1)%6]) for i in range(6)])

def cutWS(points, b1, b2, scale=1., tolerance=1e-13):
    """ 
    From given points, return only those that are inside (or on the edge of)
    the Wigner-Seitz cell of a (scale*b1, scale*b2)-based lattice.
    """
    # TODO check input dimensions?
    bvs = shortestBase46(b1, b2)
    points = np.array(points)
    for b in bvs: 
        mask = (np.tensordot(points, b, axes=(-1, 0)) <=  (scale * (1+tolerance) / 2) *np.linalg.norm(b, ord=2)**2 )
        points = points[mask]
    return points

def filledWS(b1, b2, density=10, scale=1.):
    """
    TODO doc
    TODO more intelligent generation, anisotropy balancing etc.
    """
    points = generateLattice(b1,b2,maxlayer=density*scale, include_origin=True)
    points = cutWS(points/density, np.array(b1)*scale, np.array(b2)*scale)
    return points
    
    
def reciprocalBasis1(*pargs):
    a = reduceBasisSingle(*pargs)
    return np.linalg.inv(a).T

def reciprocalBasis2pi(*pargs):
    return 2*np.pi*reciprocalBasis1(*pargs)

# TODO fill it with "points from reciprocal space" instead
def filledWS2(b1,b2, density=10, scale=1.):
    b1, b2 = reduceBasisSingle(b1,b2)
    b1r, b2r = reciprocalBasis2pi(b1,b2)
    b1l = np.linalg.norm(b1, ord=2)
    b2l = np.linalg.norm(b2, ord=2)
    b1rl = np.linalg.norm(b1r, ord=2)
    b2rl = np.linalg.norm(b2r, ord=2)
    # Black magick. Think later.™ Really. FIXME
    sicher_ratio = np.maximum(b1rl/b2rl, b2rl/b1rl) * np.maximum(b1l/b2l, b2l/b1l) # This really has to be adjusted
    points = generateLattice(b1r,b2r,maxlayer=density*scale*sicher_ratio, include_origin=True)
    points = cutWS(points*b1l/b1rl/density, b1*scale, b2*scale)
    return points


def change_basis(srcbasis, destbasis, srccoords, srccoordsaxis=-1, lattice=True):
    srcbasis = np.array(srcbasis)
    destbasis = np.array(destbasis)
    trmatrix = np.dot(np.linalg.inv(np.transpose(destbasis)), np.transpose(srcbasis))
    if lattice: # if srcbasis and destbasis are two bases of the same lattice, its elements are ints
        otrmatrix = trmatrix
        trmatrix = np.round(trmatrix)
        if not np.all(np.isclose(trmatrix, otrmatrix)):
            raise ValueError("Given srcbasis and destbasis are not bases" 
                "of the same lattice", srcbasis, destbasis, trmatrix-otrmatrix)
    destcoords = np.tensordot(srccoords, trmatrix, axes=(srccoordsaxis, -1))
    return destcoords

"""
TODO
====

- DOC!!!!!
- (nehoří) výhledově pořešit problém „hodně anisotropních“ mřížek (tj. kompensovat
rozdílné délky základních vektorů).

"""
