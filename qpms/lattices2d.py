import numpy as np
from enum import Enum

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
    TODO doc
    inputs and outputs are (2,)-shaped numpy arrays
    """
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
    return (b1,b2)

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
    # avoid obtuse angle between b1 and b2
    if b3s - b2s - b1s < eps:
        b2 = b2 + b1
        b2s = np.sum(b2 ** 2)
        b3 = b2 - b1
        b3s = np.sum(b3 ** 2)
        # This will, however, probably not happen due to the basis reduction
        print (sys.stderr, "it happened, obtuse angle!")
    if abs(b2s - b1s) < eps: # isoscele
        if abs(b3s - b1s) < eps:
            return LatticeType.EQUILATERAL_TRIANGULAR
        elif abs(b3s - 2 * b1s) < eps:
            return LatticeType.SQUARE
        else:
            return LatticeType.RHOMBIC
    elif abs(b3s - b2s - b1s) < eps:
        return LatticeType.SQUARE
    else:
        return LatticeType.OBLIQUE

def range2D(maxN, mini=1, minj=0):
    """
    "Triangle indices"
    Generates pairs of non-negative integer indices (i, j) such that
    i + j ≤ maxN, i ≥ mini, j ≥ minj.
    TODO doc and possibly different orderings
    """
    for maxn in range(min(mini, minj), maxN+1): # i + j == maxn
        for i in range(mini, maxn + 1):
            yield (i, maxn - i)

def cellWignerSeitz(b1, b2,):
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
        return np.array([solveWS(bvs[i], bvs[(i+1)%6]] for i in range(6)])





