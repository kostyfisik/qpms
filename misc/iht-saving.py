import qpms
import numpy as np
from numpy import newaxis as nx
import math
import cmath
import os
from scipy.constants import c, e as eV, hbar
s3 = math.sqrt(3)

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("omega")
#parser.add_argument("maxlayer")
args = parser.parse_args()
omega_eV = float(args.omega)

print(omega_eV)

epsilon_b = 2.3104
hexside = 375e-9
lMax = 3
maxlayer = 222
my, ny = qpms.get_mn_y(lMax)
nelem = len(my)

omega = omega_eV * eV / hbar

k_0 = omega * math.sqrt(epsilon_b) / c

output_prefix = '/tmp/diracpoints-newdata2/%d/' % maxlayer

os.makedirs(output_prefix, exist_ok=True)
qpms.hexlattice_precalc_AB_save3(file=output_prefix+str(omega_eV), lMax=lMax, k_hexside=k_0*hexside,
        maxlayer=maxlayer, savepointinfo=True)
