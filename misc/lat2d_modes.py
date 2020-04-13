#!/usr/bin/env python3

import math
from qpms.argproc import ArgParser, sfloat

ap = ArgParser(['const_real_background', 'lattice2d', 'multi_particle']) # TODO general analytical background
    
ap.add_argument("-k", nargs=2, type=sfloat, required=True, help='k vector', metavar=('K_X', 'K_Y'))
ap.add_argument("--kpi", action='store_true', help="Indicates that the k vector is given in natural units instead of SI, i.e. the arguments given by -k shall be automatically multiplied by pi / period (given by -p argument)")
ap.add_argument("--rank-tol", type=float, required=False)
ap.add_argument("-o", "--output", type=str, required=False, help='output path (if not provided, will be generated automatically)')
ap.add_argument("-t", "--rank-tolerance", type=float, default=1e11)
ap.add_argument("-c", "--min-candidates", type=int, default=1, help='always try at least this many eigenvalue candidates, even if their SVs in the rank tests are lower than rank_tolerance')
ap.add_argument("-T", "--residual-tolerance", type=float, default=2.)
ap.add_argument("-N", type=int, default="150", help="Integration contour discretisation size")


#TODO alternative specification of the contour by center and half-axes
dospec = ap.add_argument_group("Eigenvalue search area by diffracted order specification", "Specification of eigenvalue search area by diffracted order number (requires constant real refractive index for background): the integration contour 'touches' the empty lattice band specified by -b, and its axis lying on the real axis reaches '-f'-way to the next diffractive order")
dospec.add_argument("-d", "--band-index", type=int, help="Argument's absolute value determines the empty lattice band order (counted from 1), -/+ determines whether the eigenvalues are searched below/above that empty lattice band.", required=True)
dospec.add_argument("-f", "--interval-factor", type=float, default=0.1, help="Relative length of the integration ellipse axis w.r.t. the interval between two empty lattice bands; this should be be less than 1.") #TODO check
dospec.add_argument("-i", "--imaginary-aspect-ratio", type=float, default=1., help="Aspect ratio of the integration ellipse (Im/Re); this should not exceed 1/interval_factor.")

ap.add_argument("-P", "--plot", action='store_true', help="if -p not given, plot to a default path")
ap.add_argument("-O", "--plot-out", type=str, required=False, help="path to plot output (optional)")

a = ap.parse_args()

a1 = ap.direct_basis[0]
a2 = ap.direct_basis[1]

particlestr = "modes" # TODO particle string specifier or some hash, do this in argproc.py
defaultprefix = "%s_basis%gnm_%gnm__%gnm_%gnm_n%g_b%+d_k(%g_%g)um-1_cn%d" % (
            particlestr, a1[0]*1e9, a1[1]*1e9, a2[0]*1e9, a2[1]*1e9, a.refractive_index, a.band_index, a.k[0]*1e-6, a.k[1]*1e-6, a.N)


import logging
logging.basicConfig(format='%(asctime)s %(message)s', level=logging.INFO)

import numpy as np
import qpms
from qpms.cybspec import BaseSpec
from qpms.cytmatrices import CTMatrix
from qpms.qpms_c import Particle, ScatteringSystem, empty_lattice_modes_xy
from qpms.cymaterials import EpsMu, EpsMuGenerator, LorentzDrudeModel, lorentz_drude
from qpms.constants import eV, hbar
eh = eV/hbar

def inside_ellipse(point_xy, centre_xy, halfaxes_xy):
    x = point_xy[0] - centre_xy[0]
    y = point_xy[1] - centre_xy[1]
    ax = halfaxes_xy[0]
    ay = halfaxes_xy[1]
    return ((x/ax)**2 + (y/ay)**2) <= 1

beta = np.array(a.k)

if True: # TODO alternative specification of the contour by center and half-axes
    empty_freqs = empty_lattice_modes_xy(ap.background_epsmu, ap.reciprocal_basis2pi, np.array([0,0]), 1)
    empty_freqs = empty_lattice_modes_xy(ap.background_epsmu, ap.reciprocal_basis2pi, beta, (1+abs(a.band_index)) * empty_freqs[1])

    # make the frequencies in the list unique
    empty_freqs = list(empty_freqs)
    i = 0
    while i < len(empty_freqs) - 1:
        if math.isclose(empty_freqs[i], empty_freqs[i+1], rel_tol=1e-13):
            del empty_freqs[i+1]
        else:
            i += 1

    logging.info("Empty freqs: %s", str(empty_freqs))
    logging.info("Empty freqs (eV): %s", str([ff / eh for ff in empty_freqs]))
    if a.band_index > 0:
        top = empty_freqs[a.band_index]
        bottom = empty_freqs[a.band_index - 1]
        lebeta_om = bottom
    else: # a.band_index < 0
        top = empty_freqs[abs(a.band_index) - 1]
        bottom = empty_freqs[abs(a.band_index) - 2] if abs(a.band_index) > 1 else 0.
        lebeta_om = top
    #print(top,bottom,lebeta_om)
    freqradius = .5 * (top - bottom) * a.interval_factor

    centfreq = bottom + freqradius if a.band_index > 0 else top - freqradius
    if freqradius == 0:
            raise ValueError("Integration contour radius is set to zero. Are you trying to look below the lowest empty lattice band at the gamma point?")

    freqradius *= (1-1e-13) # to not totally touch the singularities

logging.info("Direct lattice basis: %s" % str(ap.direct_basis))
logging.info("Reciprocal lattice basis: %s" % str(ap.reciprocal_basis2pi))

ss, ssw = ScatteringSystem.create(ap.get_particles(), ap.background_emg, centfreq, latticebasis=ap.direct_basis)

logging.info("Finding eigenvalues around %s (= %s eV)" % (str(centfreq), str(centfreq/eh)))
logging.info("Real half-axis %s (= %s eV)" % (str(freqradius), str(freqradius/eh)))
logging.info("Imaginary half-axis %s (= %s eV)" % (str(freqradius*a.imaginary_aspect_ratio), str(freqradius*a.imaginary_aspect_ratio/eh)))

with qpms.pgsl_ignore_error(15):
    res = ss.find_modes(centfreq, freqradius, freqradius * a.imaginary_aspect_ratio,
            blochvector = a.k, contour_points = a.N, rank_tol = a.rank_tolerance,
            res_tol = a.residual_tolerance, rank_min_sel = a.min_candidates)

logging.info("Eigenfrequencies found: %s" % str(res['eigval']))
logging.info("Eigenfrequencies found (eV): %s" % str(res['eigval'] / eh))

res['inside_contour'] = inside_ellipse((res['eigval'].real, res['eigval'].imag),
         (centfreq.real, centfreq.imag), (freqradius, freqradius * a.imaginary_aspect_ratio))

#res['refractive_index_internal'] = [emg(om).n for om in res['eigval']]

#del res['omega']  If contour points are not needed...
#del res['ImTW'] # not if dbg=false anyway
outfile = defaultprefix + ".npz" if a.output is None else a.output
np.savez(outfile, meta=vars(a), empty_freqs=np.array(empty_freqs), **res)
logging.info("Saved to %s" % outfile)


if a.plot or (a.plot_out is not None):
    if len(res['eigval']) == 0:
        logging.info("No eigenvalues found; nothing to plot")
        exit(1)
    imcut = np.linspace(0, -freqradius)
    recut1 = np.sqrt(lebeta_om**2+imcut**2) # incomplete Gamma-related cut
    recut2 = np.sqrt((lebeta_om/2)**2-imcut**2) + lebeta_om/2 # odd-power-lilgamma-related cut

    import matplotlib
    matplotlib.use('pdf')
    from matplotlib import pyplot as plt

    fig = plt.figure()
    ax = fig.add_subplot(111,)
    #ax.plot(res['omega'].real/eh, res['omega'].imag/eh*1e3, ':') #res['omega'] not implemented in ScatteringSystem
    ax.add_artist(matplotlib.patches.Ellipse((centfreq.real/eh, centfreq.imag/eh*1e3),
        2*freqradius/eh, 2*freqradius*a.imaginary_aspect_ratio/eh*1e3, fill=False,
        ls=':'))
    ax.scatter(x=res['eigval'].real/eh, y=res['eigval'].imag/eh*1e3 , c = res['inside_contour']
                      )
    ax.plot(recut1/eh, imcut/eh*1e3)
    ax.plot(recut2/eh, imcut/eh*1e3)
    for i,om in enumerate(res['eigval']):
            ax.annotate(str(i), (om.real/eh, om.imag/eh*1e3))
    xmin = np.amin(res['eigval'].real)/eh
    xmax = np.amax(res['eigval'].real)/eh
    xspan = xmax-xmin
    ymin = np.amin(res['eigval'].imag)/eh*1e3
    ymax = np.amax(res['eigval'].imag)/eh*1e3
    yspan = ymax-ymin
    ax.set_xlim([xmin-.1*xspan, xmax+.1*xspan])
    ax.set_ylim([ymin-.1*yspan, ymax+.1*yspan])
    ax.set_xlabel('$\hbar \Re \omega / \mathrm{eV}$')
    ax.set_ylabel('$\hbar \Im \omega / \mathrm{meV}$')
    plotfile = defaultprefix + ".pdf" if a.plot_out is None else a.plot_out
    fig.savefig(plotfile)

exit(0)

