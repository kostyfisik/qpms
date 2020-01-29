#!/usr/bin/env python3

import math
from qpms.argproc import ArgParser


ap = ArgParser(['rectlattice2d_finite', 'single_particle', 'single_lMax', ])

ap.add_argument("-t", "--rank-tolerance", type=float, default=1e11)
ap.add_argument("-c", "--min-candidates", type=int, default=1, help='always try at least this many eigenvalue candidates, even if their SVs in the rank tests are lower than rank_tolerance')
ap.add_argument("-T", "--residual-tolerance", type=float, default=0.)

ap.add_argument("-o", "--output", type=str, required=False, help='output path (if not provided, will be generated automatically)')
#ap.add_argument("-O", "--plot-out", type=str, required=False, help="path to plot output (optional)")
#ap.add_argument("-P", "--plot", action='store_true', help="if -p not given, plot to a default path")
#ap.add_argument("-g", "--save-gradually", action='store_true', help="saves the partial result after computing each irrep")

ap.add_argument("-f", "--centre", type=complex, required=True, help='Contour centre in eV')
ap.add_argument("--ai", type=float, default=0.05, help="Contour imaginary half-axis in eV")
ap.add_argument("--ar", type=float, default=0.05, help="Contour real half-axis in eV")
ap.add_argument("-N", type=int, default="150", help="Integration contour discretisation size")

ap.add_argument("--irrep", type=str, nargs=1, default="none", help="Irrep subspace (irrep index from 0 to 7, irrep label, or 'none' for no irrep decomposition")

a=ap.parse_args()

import logging
logging.basicConfig(format='%(asctime)s %(message)s', level=logging.INFO)

Nx, Ny = a.size
px, py = a.period

particlestr = ("sph" if a.height is None else "cyl") + ("_r%gnm" % (a.radius*1e9))
if a.height is not None: particlestr += "_h%gnm" % (a.height * 1e9)
defaultprefix = "%s_p%gnmx%gnm_%dx%d_m%s_n%g_L%d_c(%s±%g±%gj)eV_cn%d" % (
    particlestr, px*1e9, py*1e9, Nx, Ny, str(a.material), a.refractive_index, a.lMax,
    str(a.centre), a.ai, a.ar, a.N)
logging.info("Default file prefix: %s" % defaultprefix)


import numpy as np
import qpms
from qpms.cybspec import BaseSpec
from qpms.cytmatrices import CTMatrix, TMatrixGenerator
from qpms.qpms_c import Particle
from qpms.cymaterials import EpsMu, EpsMuGenerator, LorentzDrudeModel, lorentz_drude
from qpms.cycommon import DebugFlags, dbgmsg_enable
from qpms import FinitePointGroup, ScatteringSystem, BesselType, eV, hbar
from qpms.symmetries import point_group_info
eh = eV/hbar

dbgmsg_enable(DebugFlags.INTEGRATION)

#Particle positions
orig_x = (np.arange(Nx/2) + (0 if (Nx % 2) else .5)) * px
orig_y = (np.arange(Ny/2) + (0 if (Ny % 2) else .5)) * py

orig_xy = np.stack(np.meshgrid(orig_x, orig_y), axis = -1)


bspec = BaseSpec(lMax = a.lMax)
medium = EpsMuGenerator(ap.background_epsmu)
particles= [Particle(orig_xy[i], ap.tmgen, bspec) for i in np.ndindex(orig_xy.shape[:-1])]

sym = FinitePointGroup(point_group_info['D2h'])
logging.info("Creating scattering system object")
ss, ssw = ScatteringSystem.create(particles, medium, sym, a.centre * eh)

if a.irrep == 'none':
    iri = None # no irrep decomposition
    irname = 'full'
    logging.info("Not performing irrep decomposition and working with the full space of dimension %d." % ss.fecv_size)
else:
    try:
        iri = int(a.irrep)
    except ValueError:
        iri = ss.irrep_names.index(a.irrep)
    irname = ss.irrep_names(iri)
    logging.info("Using irrep subspace %s (iri = %d) of dimension %d." % (irname, iri, ss.saecv_sizes[iri]))

outfile_tmp = defaultprefix + ".tmp" if a.output is None else a.output + ".tmp"

logging.info("Starting Beyn's algorithm")
results = ss.find_modes(iri=iri, omega_centre = a.centre*eh, omega_rr=a.ar*eh, omega_ri=a.ai*eh, contour_points=a.N,
        rank_tol=a.rank_tolerance, rank_min_sel=a.min_candidates, res_tol=a.residual_tolerance)
results['inside_contour'] = inside_ellipse((results['eigval'].real, results['eigval'].imag)
        (a.centre.real*eh, a.centre.imag*eh), (a.ar*eh, a.ai*eh))
results['refractive_index_internal'] = [medium(om).n for om in results['eigval']]

outfile = defaultprefix + (('_ir%s_%s.npz' % (str(iri), irname)) if iri is not None else '.npz') if a.output is None else a.output
np.savez(outfile, meta=vars(a), **results)
logging.info("Saved to %s" % outfile)

exit(0)

# TODO plots.
if a.plot or (a.plot_out is not None):
    import matplotlib
    matplotlib.use('pdf')
    from matplotlib import pyplot as plt

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(sinalpha_list, σ_ext*1e12,label='$\sigma_\mathrm{ext}$')
    ax.plot(sinalpha_list, σ_scat*1e12, label='$\sigma_\mathrm{scat}$')
    ax.plot(sinalpha_list, σ_abs*1e12, label='$\sigma_\mathrm{abs}$')
    ax.legend()
    ax.set_xlabel('$\sin\\alpha$')
    ax.set_ylabel('$\sigma/\mathrm{\mu m^2}$')
    
    plotfile = defaultprefix + ".pdf" if a.plot_out is None else a.plot_out
    fig.savefig(plotfile)

exit(0)

