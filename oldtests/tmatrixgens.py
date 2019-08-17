from qpms import *
lMax = 3
spec = BaseSpec(lMax=lMax)

ω = 1.5*eV/ℏ

gensphere_arc = TMatrixGenerator.sphere_asarc(outside=EpsMu(2.3104,1), inside=lorentz_drude['Au'], r=50e-9)
np.diag(gensphere_arc(spec,ω)[...])

