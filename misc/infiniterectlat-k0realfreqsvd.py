#!/usr/bin/env python3

import math
from qpms.argproc import ArgParser

ap = ArgParser(['rectlattice2d', 'single_particle', 'single_lMax', 'omega_seq'])
ap.add_argument("-o", "--output", type=str, required=False, help='output path (if not provided, will be generated automatically)')
ap.add_argument("-O", "--plot-out", type=str, required=False, help="path to plot output (optional)")
ap.add_argument("-P", "--plot", action='store_true', help="if -p not given, plot to a default path")
ap.add_argument("-s", "--singular_values", type=int, default=10, help="Number of singular values to plot")

a=ap.parse_args()

import logging
logging.basicConfig(format='%(asctime)s %(message)s', level=logging.INFO)

px, py = a.period

#Important! The particles are supposed to be of D2h/D4h symmetry
thegroup = 'D4h' if px == py else 'D2h'

particlestr = ("sph" if a.height is None else "cyl") + ("_r%gnm" % (a.radius*1e9))
if a.height is not None: particlestr += "_h%gnm" % (a.height * 1e9)
defaultprefix = "%s_p%gnmx%gnm_m%s_n%g_f(%g..%g..%g)eV_L%d_SVGamma" % (
    particlestr, px*1e9, py*1e9, str(a.material), a.refractive_index, *(a.eV_seq), a.lMax)
logging.info("Default file prefix: %s" % defaultprefix)


import numpy as np
import qpms
from qpms.cybspec import BaseSpec
from qpms.cytmatrices import CTMatrix, TMatrixGenerator
from qpms.qpms_c import Particle, pgsl_ignore_error
from qpms.cymaterials import EpsMu, EpsMuGenerator, LorentzDrudeModel, lorentz_drude
from qpms.cycommon import DebugFlags, dbgmsg_enable
from qpms import FinitePointGroup, ScatteringSystem, BesselType, eV, hbar
from qpms.cyunitcell import unitcell
#from qpms.symmetries import point_group_info # TODO
eh = eV/hbar

# not used; TODO:
irrep_labels = {"B2''":"$B_2''$",
                "B2'":"$B_2'$",
                "A1''":"$A_1''$",
                "A1'":"$A_1'$",
                "A2''":"$A_2''$",
                "B1''":"$B_1''$",
                "A2'":"$A_2'$", 
                "B1'":"$B_1'$"}

dbgmsg_enable(DebugFlags.INTEGRATION)

a1 = ap.direct_basis[0]
a2 = ap.direct_basis[1]

#Particle positions
orig_x = [0]
orig_y = [0]
orig_xy = np.stack(np.meshgrid(orig_x,orig_y),axis=-1)

omegas = ap.omegas

logging.info("%d frequencies from %g to %g eV" % (len(omegas), omegas[0]/eh, omegas[-1]/eh))

bspec = BaseSpec(lMax = a.lMax)
nelem = len(bspec)
# The parameters here should probably be changed (needs a better qpms_c.Particle implementation)
pp = Particle(orig_xy[0][0], tmgen=ap.tmgen, bspec=bspec)
par = [pp]

u = unitcell(a1, a2, par, refractive_index=a.refractive_index)
eta = (np.pi / u.Area)**.5

wavenumbers = np.empty(omegas.shape)
SVs = np.empty(omegas.shape+(nelem,))
beta = np.array([0.,0.])
for i, omega in enumerate(omegas):
    wavenumbers[i] = ap.background_epsmu.k(omega).real # Currently, ScatteringSystem does not "remember" frequency nor wavenumber
    with pgsl_ignore_error(15): # avoid gsl crashing on underflow; maybe not needed
        ImTW = u.evaluate_ImTW(eta, omega, beta)
    SVs[i] = np.linalg.svd(ImTW, compute_uv = False)

outfile = defaultprefix + ".npz" if a.output is None else a.output
np.savez(outfile, meta=vars(a), omegas=omegas, wavenumbers=wavenumbers, SVs=SVs, unitcell_area=u.Area
       )
logging.info("Saved to %s" % outfile)


if a.plot or (a.plot_out is not None):
    import matplotlib
    matplotlib.use('pdf')
    from matplotlib import pyplot as plt

    fig = plt.figure()
    ax = fig.add_subplot(111)
    for i in range(a.singular_values):
        ax.plot(omegas/eh, SVs[:,-1-i])
    ax.set_xlabel('$\hbar \omega / \mathrm{eV}$')
    ax.set_ylabel('Singular values')
    
    plotfile = defaultprefix + ".pdf" if a.plot_out is None else a.plot_out
    fig.savefig(plotfile)

exit(0)

