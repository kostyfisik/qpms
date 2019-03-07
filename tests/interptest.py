
from qpms import Particle, CTMatrix, BaseSpec, FinitePointGroup, ScatteringSystem, TMatrixInterpolator
bspec1 = BaseSpec(lMax = 1)

tmfile = '/home/mmn/tmp/cylinder_50nm_lMax4_cleaned.TMatrix'
interp = TMatrixInterpolator(tmfile, bspec1)
interp(3)
