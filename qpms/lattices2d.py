import numpy as np
import warnings
from enum import Enum

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
    return(b1,b2)

def orderedReducedBasis(b1, b2):
    ''' blah blab blah
    |b1| is still the shortest possible basis vector,
        but if there would be obtuse angle between b1 and b2, b2 - b1 is returned
        in place of the original b2. In other words, b1, b2 and b2-b1 are 
    '''
    b1, b2 = reduceBasisSingle(b1,b2)
    
    if b3s - b2s - b1s > eps: # obtuse angle between b1 and b2
        pass
    pass
    #-------- zde jsem skončil ------------


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
        warnings.warn("obtuse angle between reduced basis vectors, the lattice type identification might is not well tested.")
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
    for maxn in range(min(mini, minj, minN), maxN+1): # i + j == maxn
        for i in range(mini, maxn + 1):
            yield (i, maxn - i)

            
def generateLattice(b1, b2, maxlayer=5, include_origin=False, order='leaves'):
    b1, b2 = reduceBasisSingle(b1, b2)
    latticeType = classifyLatticeSingle(b1, b2)


    if latticeType is LatticeType.RECTANGULAR or latticeType is LatticeType.SQUARE:
        bvs = (b1, b2, -b1, -b2)
    else:
            # Avoid obtuse angle between b1 and b2. TODO This should be yet thoroughly tested.
        if is_obtuse(b1,b2):
            b3 = b2
            b2 = b2 + b1
            # N. B. now the assumption |b3| >= |b2| is no longer valid
            warnings.warn("obtuse angle between reduced basis vectors, the lattice generation might is not well tested.")
        else:
            b3 = b2 - b1
        bvs = (b1, b2, b3, -b1, -b2, -b3)
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
    latticeType = classifyLaticeSingle(b1, b2)
    if latticeType is LatticeType.RECTANGULAR or latticeType is LatticeType.SQUARE:
        return np.array( (
            (+b1+b2),
            (+b2-b1),
            (-b1-b2),
            (-b2+b1),
        )) / 2
    else:
        b3 = b2 - b1
        bvs = (b1, b2, b3, -b1, -b2, -b3)
        return np.array([solveWS(bvs[i], bvs[(i+1)%6]) for i in range(6)])

def cutWS(points, b1, b2, scale=1.):
    """ 
    From given points, return only those that are inside (or on the edge of)
    the Wigner-Seitz cell of a (scale*b1, scale*b2)-based lattice.
    """
    # TODO check input dimensions?
    b1, b2 = reduceBasisSingle(b1, b2)
    b3 = b2 - b1
    bvs = (b1, b2, b3, -b1, -b2, -b3)
    points = np.array(points)
    for b in bvs: 
        mask = (np.tensordot(points, b, axes=(-1, 0)) <= np.linalg.norm(b, ord=2) * scale/2)
        points = points[mask]
    return points

def filledWS(b1, b2, density=10, scale=1.):
    """
    TODO doc
    TODO more intelligent generation, anisotropy balancing etc.
    """
    b1, b2 = reduceBasisSingle(b1, b2)
    pass
    
    

def reciprocalBasis(a1, a2):
    pass
