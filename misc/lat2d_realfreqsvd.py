#!/usr/bin/env python3

import math
from qpms.argproc import ArgParser, sfloat

ap = ArgParser(['background', 'lattice2d', 'multi_particle',  'omega_seq'])

ap.add_argument("-k", nargs=2, type=sfloat, required=True, help='k vector', metavar=('K_X', 'K_Y'))
ap.add_argument("--kpi", action='store_true', help="Indicates that the k vector is given in natural units instead of SI, i.e. the arguments given by -k shall be automatically multiplied by pi / period (given by -p argument)")

ap.add_argument("-g", "--little-group", type=str, default="trivial_g", help="Little group for subspace irrep classification", action="store")

ap.add_argument("-o", "--output", type=str, required=False, help='output path (if not provided, will be generated automatically)')
ap.add_argument("-O", "--plot-out", type=str, required=False, help="path to plot output (optional)")
ap.add_argument("-P", "--plot", action='store_true', help="if -p not given, plot to a default path")
ap.add_argument("-s", "--singular_values", type=int, default=10, help="Number of singular values to plot")

a=ap.parse_args()

import logging
logging.basicConfig(format='%(asctime)s %(message)s', level=logging.INFO)


#Important! The particles are supposed to be of D2h/D4h symmetry
# thegroup = 'D4h' if px == py and not a.D2 else 'D2h'

a1 = ap.direct_basis[0]
a2 = ap.direct_basis[1]

particlestr = "svdinterval" # TODO particle string specifier or some hash, do this in argproc.py
defaultprefix = "%s_basis%gnm_%gnm__%gnm_%gnm_f(%g..%g..%g)eV_k%g_%g" % (
    particlestr, a1[0]*1e9, a1[1]*1e9, a2[0]*1e9, a2[1]*1e9, *(a.eV_seq), ap.k[0], ap.k[1])
logging.info("Default file prefix: %s" % defaultprefix)


import numpy as np
import qpms
import warnings
from qpms.cybspec import BaseSpec
from qpms.cytmatrices import CTMatrix, TMatrixGenerator
from qpms.qpms_c import Particle, pgsl_ignore_error, empty_lattice_modes_xy
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



omegas = ap.omegas

logging.info("%d frequencies from %g to %g eV" % (len(omegas), omegas[0]/eh, omegas[-1]/eh))

particles = ap.get_particles()

ss, ssw = ScatteringSystem.create(particles, ap.background_emg, omegas[0], latticebasis=ap.direct_basis)
k = np.array([ap.k[0], ap.k[1], 0])
# Auxillary finite scattering system for irrep decomposition, quite a hack
ss1, ssw1 = ScatteringSystem.create(particles, ap.background_emg, omegas[0],sym=FinitePointGroup(point_group_info[ap.little_group]))

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
                label= None if i else irrep_labels.get(ss1.irrep_names[iri], ss1.irrep_names[iri]), 
                **cargs)
    ax.set_ylim([0,1.1])
    if hasattr(ap, "background_epsmu"):
        xlim = ax.get_xlim()
        omegas_empty = empty_lattice_modes_xy(ap.background_epsmu, ap.reciprocal_basis2pi, k, omegas[-1])
        for om in omegas_empty:
            if om/eh > xlim[0] and om/eh < xlim[1]:
                ax.axvline(om/eh, ls=':')
    ax.set_xlabel('$\hbar \omega / \mathrm{eV}$')
    ax.set_ylabel('Singular values')
    ax.legend()
    
    plotfile = defaultprefix + ".pdf" if a.plot_out is None else a.plot_out
    fig.savefig(plotfile)

exit(0)

