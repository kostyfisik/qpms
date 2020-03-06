#!/usr/bin/env python3

import math
from qpms.argproc import ArgParser

ap = ArgParser(['rectlattice2d', 'single_particle', 'single_lMax', 'single_omega'])
ap.add_argument("-k", '--kx-lim', nargs=2, type=float, required=True, help='k vector', metavar=('KX_MIN', 'KX_MAX'))
# ap.add_argument("--kpi", action='store_true', help="Indicates that the k vector is given in natural units instead of SI, i.e. the arguments given by -k shall be automatically multiplied by pi / period (given by -p argument)")
ap.add_argument("-o", "--output", type=str, required=False, help='output path (if not provided, will be generated automatically)')
ap.add_argument("-N", type=int, default="151", help="Number of angles")
ap.add_argument("-O", "--plot-out", type=str, required=False, help="path to plot output (optional)")
ap.add_argument("-P", "--plot", action='store_true', help="if -p not given, plot to a default path")
#ap.add_argument("-g", "--save-gradually", action='store_true', help="saves the partial result after computing each irrep")


a=ap.parse_args()

import logging
logging.basicConfig(format='%(asctime)s %(message)s', level=logging.INFO)

px, py = a.period

particlestr = ("sph" if a.height is None else "cyl") + ("_r%gnm" % (a.radius*1e9))
if a.height is not None: particlestr += "_h%gnm" % (a.height * 1e9)
defaultprefix = "%s_p%gnmx%gnm_m%s_n%g_angles(%g_%g)_Ey_f%geV_L%d_cn%d" % (
    particlestr, px*1e9, py*1e9, str(a.material), a.refractive_index, a.kx_lim[0], a.kx_lim[1], a.eV, a.lMax, a.N)
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
eh = eV/hbar

dbgmsg_enable(DebugFlags.INTEGRATION)

a1 = ap.direct_basis[0]
a2 = ap.direct_basis[1]

#Particle positions
orig_x = [0]
orig_y = [0]
orig_xy = np.stack(np.meshgrid(orig_x,orig_y),axis=-1)

omega = ap.omega

bspec = BaseSpec(lMax = a.lMax)
# The parameters here should probably be changed (needs a better qpms_c.Particle implementation)
pp = Particle(orig_xy[0][0], ap.tmgen, bspec=bspec)
par = [pp]


ss, ssw = ScatteringSystem.create(par, ap.background_emg, omega, latticebasis = ap.direct_basis)

if ssw.wavenumber.imag != 0:
    warnings.warn("The background medium wavenumber has non-zero imaginary part. Don't expect meaningful results for cross sections.")
wavenumber = ssw.wavenumber.real 

sinalpha_list = np.linspace(a.kx_lim[0],a.kx_lim[1],a.N)

# Plane wave data
E_cart_list = np.empty((a.N,3), dtype=complex)
E_cart_list[:,:] = np.array((0,1,0))[None,:]
k_cart_list = np.empty((a.N,3), dtype=float)
k_cart_list[:,0] = sinalpha_list
k_cart_list[:,1] = 0
k_cart_list[:,2] = np.sqrt(1-sinalpha_list**2)
k_cart_list *= wavenumber

σ_ext_list = np.empty((a.N,), dtype=float)
σ_scat_list = np.empty((a.N,), dtype=float)

with pgsl_ignore_error(15): # avoid gsl crashing on underflow
    for j in range(a.N):
        k_cart = k_cart_list[j]
        blochvector = (k_cart[0], k_cart[1], 0)
        # the following two could be calculated only once, but probably not a big deal
        LU = ssw.scatter_solver(k=blochvector)
        ã = ss.planewave_full(k_cart=k_cart, E_cart=E_cart_list[j])
        Tã = ssw.apply_Tmatrices_full(ã)
        f = LU(Tã)

        σ_ext_list[j] = -np.vdot(ã, f).real/wavenumber**2
        translation_matrix = ssw.translation_matrix_full(blochvector=blochvector) + np.eye(ss.fecv_size)
        σ_scat_list[j] = np.vdot(f,np.dot(translation_matrix, f)).real/wavenumber**2

σ_abs_list = σ_ext_list - σ_scat_list

outfile = defaultprefix + ".npz" if a.output is None else a.output
np.savez(outfile, meta=vars(a), sinalpha=sinalpha_list, k_cart = k_cart_list, E_cart=E_cart_list, σ_ext=σ_ext_list,σ_abs=σ_abs_list,σ_scat=σ_scat_list, omega=omega, wavenumber=wavenumber, unitcell_area=ss.unitcell_volume
       )
logging.info("Saved to %s" % outfile)


if a.plot or (a.plot_out is not None):
    import matplotlib
    matplotlib.use('pdf')
    from matplotlib import pyplot as plt

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(sinalpha_list, σ_ext_list*1e12,label='$\sigma_\mathrm{ext}$')
    ax.plot(sinalpha_list, σ_scat_list*1e12, label='$\sigma_\mathrm{scat}$')
    ax.plot(sinalpha_list, σ_abs_list*1e12, label='$\sigma_\mathrm{abs}$')
    ax.legend()
    ax.set_xlabel('$\sin\\alpha$')
    ax.set_ylabel('$\sigma/\mathrm{\mu m^2}$')
    
    plotfile = defaultprefix + ".pdf" if a.plot_out is None else a.plot_out
    fig.savefig(plotfile)

exit(0)

