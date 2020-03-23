#!/usr/bin/env python3
"""
Tests whether direct evaluation of VSWFs gives the same result as their
representation in terms of translation operators and regular electric dipole
waves at origin
"""

from qpms import Particle, CTMatrix, lorentz_drude, EpsMuGenerator, TMatrixGenerator, BaseSpec, FinitePointGroup, ScatteringSystem, TMatrixInterpolator, EpsMu, dbgmsg_enable, dbgmsg_disable, dbgmsg_active, BesselType,eV, hbar
from qpms.qpms_c import set_gsl_pythonic_error_handling
from qpms.symmetries import point_group_info
import numpy as np
eh = eV/hbar
np.random.seed(666)
dbgmsg_enable(2)

part_radius = 80e-9
p = 1580e-9

set_gsl_pythonic_error_handling()

sym = FinitePointGroup(point_group_info['D4h'])
bspec1 = BaseSpec(lMax=3)
medium=EpsMuGenerator(EpsMu(1.52**2))
t1 = TMatrixGenerator.sphere(medium, EpsMuGenerator(lorentz_drude['Au']), r=part_radius)
p1 = Particle((0,0,0),t1,bspec=bspec1)
ss, ssw = ScatteringSystem.create([p1], EpsMuGenerator(EpsMu(1.52**2)), 1.4*eh, sym)


points = np.random.random((100,3)) * p
points = points[np.linalg.norm(points, axis=-1) > part_radius]
t,l,m = bspec1.tlm()

fails=0

for i in range(ss.fecv_size):
    fvc = np.zeros((ss.fecv_size,), dtype=complex)
    fvc[i] = 1
    
    E = ssw.scattered_E(fvc, points)
    E_alt = ssw.scattered_E(fvc, points,alt=True)
    diff = abs(E-E_alt)
    reldiffavg = np.nanmean(diff/(abs(E)+abs(E_alt)))
    fail = reldiffavg > 1e-3
    fails += fail
                                
    print('E' if t[i] == 2 else 'M', l[i], m[i], np.amax(diff), reldiffavg, 'FAIL!' if fail else 'OK')

exit(fails)
