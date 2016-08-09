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
lengthOrdersOfMagnitude = [2.**i for i in range(-15,10,2)]

class PlaneWaveDecompositionTests(unittest.TestCase):
    """
    covers plane_pq_y and vswf_yr1
    """
    def testRandomInc(self):
        # The "maximum" argument of the Bessel's functions, i.e. maximum wave number times the distance,
        # for the "locally strongly varying fields"
        maxx = 10
        rfailtol = 0.01 # how much of the randomized test will be tolerated
        lMax = 80 # To which order we decompose the waves
        rtol = 1e-5 # relative required precision
        atol = 1. # absolute tolerance, does not really play a role
        nsamples = 4 # (frequency, direction, polarisation) samples per order of magnitude and test
        npoints = 15 # points to evaluate per sample

    
        failcounter = 0
        passcounter = 0
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
                    #self.assertTrue(np.allclose(planewave_2_p, planewave_1[i], rtol=rtol, atol=atol))
                    if not np.allclose(planewave_2_p, planewave_1[i], rtol=rtol, atol=atol):
                        False and warnings.warn('Planewave expansion test not passed; r = '
                                +str(testpoints[i])+', k = '+str(k[s])
                                +', E_0 = '+str(E_0[s])+', (original) E = '
                                +str(planewave_1[i])+', (reexpanded) E = '
                                +str(planewave_2_p)
                                +', x = '+str(np.dot(testpoints[i],k[s]))
                                +'; distance = '
                                +str(np.linalg.norm(planewave_1[i]-planewave_2_p))
                                +', required relative precision = '
                                +str(rtol)+'.')
                        failcounter += 1
                    else:
                        passcounter += 1
        self.assertLess(failcounter / (failcounter + passcounter), rfailtol, 
                '%d / %d (%.2e) randomized numerical tests failed (tolerance %.2e)' 
                % (failcounter, failcounter + passcounter, 
                    failcounter / (failcounter + passcounter), rfailtol))
        return 

    def testCornerCases(self):
        pass


class SphericalWaveTranslationTests(unittest.TestCase):
    def testRandom1to1(self):
        # The "maximum" argument of the Bessel's functions, i.e. maximum wave number times the distance,
        # for the "locally strongly varying fields"
        maxx = 10
        rfailtol = 0.01 # how much of the randomized test fail proportion will be tolerated
        lMax = 50 # To which order we decompose the waves
        lMax_outgoing = 4 # To which order we try the outgoing waves
        rtol = 1e-5 # relative required precision
        atol = 1. # absolute tolerance, does not really play a role
        nsamples = 4 # frequency samples per order of magnitude and test
        npoints = 15 # points to evaluate per frequency and center
     
        ncentres = 3 # number of spherical coordinate centres between which the translations are to be made
        maxxd = 2000 # the center position standard deviation

        failcounter = 0
        passcounter = 0
        my, ny = qpms.get_mn_y(lMax)
        nelem_full = len(my)
        nelem_out = lMax_outgoing * (lMax_outgoing + 2)
        for oom in lengthOrdersOfMagnitude:
            centres = oom * maxxd * np.random.randn(ncentres, 3)
            ksizs = np.random.randn(nsamples)
            for ksiz in ksizs:
                for i in range(ncentres): # "source"
                    Rs = centres[i]
                    testr = oom * maxx * np.random.randn(npoints, 3)
                    for j in range(ncentres): # "destination"
                        if j == i:
                            continue
                        Rd = centres[j]
                        shift = Rd - Rs
                        shift_sph = qpms.cart2sph(shift)
                        shift_kr = ksiz * shift_sph[0]
                        shift_theta = shift_sph[1]
                        shift_phi = shift_sph[2]

                        A_yd_ys = np.empty((nelem_full,nelem_out), dtype = np.complex_)
                        B_yd_ys = np.empty((nelem_full,nelem_out), dtype = np.complex_)
                        for yd in range(nelem_full):
                            for ys in range(nelem_out):
                                A_yd_ys[yd, ys] = qpms.Ã(my[yd],ny[yd],my[ys],ny[ys],shift_kr, shift_theta, shift_phi, True, 1) 
                                B_yd_ys[yd, ys] = qpms.B̃(my[yd],ny[yd],my[ys],ny[ys],shift_kr, shift_phi, shift_phi, True, 1)
                        for r in testr:
                            sph_ssys = qpms.cart2sph(r+Rd-Rs)
                            M_ssys, N_ssys = qpms.vswf_yr1(np.array([ksiz * sph_ssys[0], sph_ssys[1], sph_ssys[2]]), lMax_outgoing, J=1)
                            sph_dsys = qpms.cart2sph(r)
                            M_dsys, N_dsys = qpms.vswf_yr1(np.array([ksiz * sph_dsys[0], sph_dsys[1], sph_dsys[2]]), lMax, J=1)
                            for ys in range(nelem_out): 
                                # Electrical waves
                                E_1 = -1j*qpms.sph_loccart2cart(N_ssys[ys], sph_ssys)
                                E_2 = -1j*qpms.sph_loccart2cart(np.dot(A_yd_ys[:,ys],N_dsys)+np.dot(B_yd_ys[:,ys],M_dsys),sph_dsys)
                                if not np.allclose(E_1, E_2, rtol=rtol, atol=atol):
                                    failcounter += 1
                                else:
                                    passcounter += 1
                                # Magnetic waves
                                E_1 = -1j*qpms.sph_loccart2cart(M_ssys[ys], sph_ssys)
                                E_2 = -1j*qpms.sph_loccart2cart(np.dot(A_yd_ys[:,ys],M_dsys)+np.dot(B_yd_ys[:,ys],N_dsys),sph_dsys)
                                if not np.allclose(E_1, E_2, rtol=rtol, atol=atol):
                                    failcounter += 1
                                else:
                                    passcounter += 1
        self.assertLess(failcounter / (failcounter + passcounter), rfailtol, 
                '%d / %d (%.2e) randomized numerical tests failed (tolerance %.2e)' 
                % (failcounter, failcounter + passcounter, 
                    failcounter / (failcounter + passcounter), rfailtol))
        return 

    def testRandom3to1(self):
        pass

def main():
    unittest.main()

if __name__ == '__main__':
    main()

