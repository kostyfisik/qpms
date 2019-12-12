# unit conversions, mostly for standalone usage
# TODO avoid importing the "heavy" qpms parts
from scipy.constants import epsilon_0 as ε_0, c, pi as π, e as eV, hbar, hbar as ℏ, mu_0 as μ_0
pi = π
μm = 1e-6
nm = 1e-9
# "SCUFF FREQUENCY UNIT"
SU = 3e14
SCUFF_OMEGAUNIT = SU

#freqs = freqs_weirdunits * c / μm

def nm2eV(lambda_nm, n = 1):
    return 2*π*c/(n * lambda_nm * nm) / (eV / ℏ)

def eV2nm(e_eV, n = 1):
    return 2*π*c/(n * e_eV * (eV/ℏ)) / nm

def SU2eV(e_SU):
    '''
    Scuff frequency units to eV.
    '''
    return e_SU * SU / (eV / ℏ)

def eV2SU(e_eV):
    return e_eV * (eV / ℏ) / SU

def SU2nm(e_SU, n = 1):
    return 2*π*c/(n * e_SU * SU) / nm

def nm2SU(lambda_nm, n = 1):
    return 2*π*c/(n * lambda_nm * nm) / SU

