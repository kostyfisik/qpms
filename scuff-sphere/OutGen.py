import numpy as np
from math import sin, cos, pi
from scipy.constants import c, e as eV, hbar as ℏ
# Scuff OmegaFile has units of c/1μm
ωlist = np.linspace(0, 5*eV/ℏ / (c/1e-6),100)
with open('OmegaList', 'w') as olf:
    for ω in ωlist:
        olf.write(str(ω) + '\n')

# Create a ring of N points where the field will be evaluated
R = 0.008 # μm
N = 30
with open('EPFile_' + str(N) + 'ring-xz', 'w') as pf:
    for fi in range(N):
        pf.write(str(R*sin(fi)) + '\t0\t' + str(R*cos(fi)) + '\n')

