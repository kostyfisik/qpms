from qpms import *
import numpy as np
R = 40e-9
lMax = 2
ω = 1.5*eV/ℏ

spec = BaseSpec(lMax = lMax)

inside = EpsMuGenerator(lorentz_drude['Au'])
outside = EpsMuGenerator(EpsMu(2.3104,1))

gensphere_arc = TMatrixGenerator.sphere_asarc(outside=outside, inside=inside, r=R, lMax_extend=lMax)

np.set_printoptions(precision=3, suppress=True,linewidth=1000)

QT = gensphere_arc.Q_transposed(ω, spec.norm)
RT = gensphere_arc.R_transposed(ω, spec.norm)
T = gensphere_arc(spec,ω)


QT_corrected = np.array(QT)
RT_corrected = np.array(RT)
QT_corrected[8:,:8] = 0
QT_corrected[:8,8:] = 0
RT_corrected[8:,:8] = 0
RT_corrected[:8,8:] = 0
T_corrected = np.dot(np.linalg.inv(QT_corrected), RT_corrected)

print("QT:")
print(QT)

print("RT:")
print(RT)

print("T:")
print(T[:])

print("T_corrected:")
print(T_corrected)



