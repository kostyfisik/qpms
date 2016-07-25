"""
Unit tests for qpms_p
=====================

Covered functions
-----------------
plane_pq_y vs. vswf_yr1

Not covered
-----------
Everything else

"""

import unittest
import qpms
import numpy as np
from numpy import newaxis as ň
import warnings

# Some constants go here.

# The "maximum" argument of the Bessel's functions, i.e. maximum wave number times the distance,
# for the "locally strongly varying fields"
maxx = 3
# The "maximum" argument of the Bessel's function for reexpansion of the waves into regular waves
# in another center
maxxd = 2000
lMax = 50 # To which order we decompose the waves

lengthOrdersOfMagnitude = [2.**i for i in range(-15,10)]
nsamples = 4 # (frequency, direction, polarisation) samples per order of magnitude and test
npoints = 40 # points to evaluate per sample
rtol = 1e-7 # relative required precision
atol = 1. # absolute tolerance, does not really play a role

class PlaneWaveDecompositionTests(unittest.TestCase):
    """
    covers plane_pq_y and vswf_yr1
    """
    def testRandomInc(self):
        for oom in lengthOrdersOfMagnitude:
            k = np.random.randn(nsamples, 3) / oom
            ksiz = np.linalg.norm(k, axis=-1)
            kdir = k / ksiz[...,ň]
            E_0 = np.cross(np.random.randn(nsamples, 3), k) * oom # ensure orthogonality
            for s in range(nsamples):
                testpoints = oom * maxx * np.random.randn(npoints, 3)
                p, q = qpms.plane_pq_y(lMax, k[s], E_0[s])
                planewave_1 =  np.exp(1j*np.dot(testpoints,k[s]))[:,ň] * E_0[s,:]
                for i in range(npoints):
                    sph = qpms.cart2sph(ksiz[s]*testpoints[i])
                    M̃_y, Ñ_y = qpms.vswf_yr1(sph, lMax, 1)
                    planewave_2_p = -1j*qpms.sph_loccart2cart(np.dot(p,Ñ_y)+np.dot(q,M̃_y),sph)
                    self.assertTrue(np.allclose(planewave_2_p, planewave_1[i], rtol=rtol, atol=atol))
#                    if not np.allclose(planewave_2_p, planewave_1[i], rtol=rtol, atol=atol):
#                        warnings.warn('Planewave expansion test not passed; r = '
#                                +str(testpoints[i])+', k = '+str(k[s])
#                                +', E_0 = '+str(E_0[s])+', (original) E = '
#                                +str(planewave_1[i])+', (reexpanded) E = '
#                                +str(planewave_2_p)
#                                +', x = '+str(np.dot(testpoints[i],k[s]))
#                                +'; distance = '
#                                +str(np.linalg.norm(planewave_1[i]-planewave_2_p))
#                                +', required relative precision = '
#                                +str(relprecision)+'.')
        return

    def testCornerCases(self):
        pass


class SphericalWaveTranslationTests(unittest.TestCase):

    def sometest(self):
        pass

def main():
    unittest.main()

if __name__ == '__main__':
    main()

