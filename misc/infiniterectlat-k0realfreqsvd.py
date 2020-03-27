#!/usr/bin/env python3

import math
from qpms.argproc import ArgParser

ap = ArgParser(['rectlattice2d', 'single_particle', 'single_lMax', 'omega_seq'])
ap.add_argument("-o", "--output", type=str, required=False, help='output path (if not provided, will be generated automatically)')
ap.add_argument("-O", "--plot-out", type=str, required=False, help="path to plot output (optional)")
ap.add_argument("-P", "--plot", action='store_true', help="if -p not given, plot to a default path")
ap.add_argument("-s", "--singular_values", type=int, default=10, help="Number of singular values to plot")
ap.add_argument("--D2", action='store_true', help="Use D2h symmetry even if the x and y periods are equal")

a=ap.parse_args()

import logging
logging.basicConfig(format='%(asctime)s %(message)s', level=logging.INFO)

px, py = a.period

#Important! The particles are supposed to be of D2h/D4h symmetry
thegroup = 'D4h' if px == py and not a.D2 else 'D2h'

particlestr = ("sph" if a.height is None else "cyl") + ("_r%gnm" % (a.radius*1e9))
if a.height is not None: particlestr += "_h%gnm" % (a.height * 1e9)
defaultprefix = "%s_p%gnmx%gnm_m%s_n%g_f(%g..%g..%g)eV_L%d_SVGamma" % (
    particlestr, px*1e9, py*1e9, str(a.material), a.refractive_index, *(a.eV_seq), a.lMax)
logging.info("Default file prefix: %s" % defaultprefix)


import numpy as np
import qpms
import warnings
from qpms.cybspec import BaseSpec
from qpms.cytmatrices import CTMatrix, TMatrixGenerator
from qpms.qpms_c import Particle, pgsl_ignore_error
from qpms.cymaterials import EpsMu, EpsMuGenerator, LorentzDrudeModel, lorentz_drude
from qpms.cycommon import DebugFlags, dbgmsg_enable
from qpms import FinitePointGroup, ScatteringSystem, BesselType, eV, hbar
from qpms.symmetries import point_group_info 
eh = eV/hbar

# not used; TODO:
irrep_labels = {"B2''":"$B_2''$",
                "B2'":"$B_2'$",
                "A1''":"$A_1''$",
                "A1'":"$A_1'$",
                "A2''":"$A_2''$",
                "B1''":"$B_1''$",
                "A2'":"$A_2'$", 
                "B1'":"$B_1'$",
                "E'":"$E'$",
                "E''":"$E''$",}

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
pp = Particle(orig_xy[0][0], ap.tmgen, bspec=bspec)

ss, ssw = ScatteringSystem.create([pp], ap.background_emg, omegas[0], latticebasis=ap.direct_basis)
k = np.array([0.,0.,0])
# Auxillary finite scattering system for irrep decomposition, quite a hack
ss1, ssw1 = ScatteringSystem.create([pp], ap.background_emg, omegas[0],sym=FinitePointGroup(point_group_info[thegroup]))

wavenumbers = np.empty(omegas.shape)
SVs = [None] * ss1.nirreps
for iri in range(ss1.nirreps):
    SVs[iri] = np.empty(omegas.shape+(ss1.saecv_sizes[iri],))
for i, omega in enumerate(omegas):
    ssw = ss(omega)
    wavenumbers[i] = ssw.wavenumber.real
    if ssw.wavenumber.imag: 
        warnings.warn("Non-zero imaginary wavenumber encountered")
    with pgsl_ignore_error(15): # avoid gsl crashing on underflow; maybe not needed
        ImTW = ssw.modeproblem_matrix_full(k)
    for iri in range(ss1.nirreps):
        ImTW_packed = ss1.pack_matrix(ImTW, iri)
        SVs[iri][i] = np.linalg.svd(ImTW_packed, compute_uv = False)

outfile = defaultprefix + ".npz" if a.output is None else a.output
np.savez(outfile, meta=vars(a), omegas=omegas, wavenumbers=wavenumbers, SVs=np.concatenate(SVs, axis=-1), irrep_names=ss1.irrep_names, irrep_sizes=ss1.saecv_sizes, unitcell_area=ss.unitcell_volume
       )
logging.info("Saved to %s" % outfile)


if a.plot or (a.plot_out is not None):
    import matplotlib
    matplotlib.use('pdf')
    from matplotlib import pyplot as plt

    fig = plt.figure()
    ax = fig.add_subplot(111)
    cc = plt.rcParams['axes.prop_cycle']()
    for iri in range(ss1.nirreps):
        cargs = next(cc)
        nlines = min(a.singular_values, ss1.saecv_sizes[iri])
        for i in range(nlines):
            ax.plot(omegas/eh, SVs[iri][:,-1-i],
                label= None if i else irrep_labels[ss1.irrep_names[iri]], 
                **cargs)
    ax.set_ylim([0,1.1])
    ax.set_xlabel('$\hbar \omega / \mathrm{eV}$')
    ax.set_ylabel('Singular values')
    ax.legend()
    
    plotfile = defaultprefix + ".pdf" if a.plot_out is None else a.plot_out
    fig.savefig(plotfile)

exit(0)

