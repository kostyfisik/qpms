#!/usr/bin/env python3

import math
from qpms.argproc import ArgParser


ap = ArgParser(['single_particle', 'single_omega', 'single_lMax'])
ap.add_argument("-p", "--period", type=float, required=True, help='square lattice period')
ap.add_argument("--Nx", type=int, required=True, help='Array size x')
ap.add_argument("--Ny", type=int, required=True, help='Array size y')
ap.add_argument("-k", '--kx-lim', nargs=2, type=float, required=True, help='k vector', metavar=('KX_MIN', 'KX_MAX'))
# ap.add_argument("--kpi", action='store_true', help="Indicates that the k vector is given in natural units instead of SI, i.e. the arguments given by -k shall be automatically multiplied by pi / period (given by -p argument)")
ap.add_argument("--rank-tol", type=float, required=False)
ap.add_argument("-o", "--output", type=str, required=False, help='output path (if not provided, will be generated automatically)')
ap.add_argument("-N", type=int, default="151", help="Number of angles")
ap.add_argument("-O", "--plot-out", type=str, required=False, help="path to plot output (optional)")
ap.add_argument("-P", "--plot", action='store_true', help="if -p not given, plot to a default path")
ap.add_argument("-g", "--save-gradually", action='store_true', help="saves the partial result after computing each irrep")


a=ap.parse_args()

import logging
logging.basicConfig(format='%(asctime)s %(message)s', level=logging.INFO)


particlestr = ("sph" if a.height is None else "cyl") + ("_r%gnm" % (a.radius*1e9))
if a.height is not None: particlestr += "_h%gnm" % (a.height * 1e9)
defaultprefix = "%s_p%gnm_%dx%d_m%s_n%g_angles(%g_%g)_Ey_f%geV_L%d_cn%d" % (
    particlestr, a.period*1e9, a.Nx, a.Ny, str(a.material), a.refractive_index, a.kx_lim[0], a.kx_lim[1], a.eV, a.lMax, a.N)
logging.info("Dafault file prefix: %s" % defaultprefix)


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

px=a.period
py=a.period

#Particle positions
orig_x = (np.arange(a.Nx/2) + (0 if (a.Nx % 2) else .5)) * px
orig_y = (np.arange(a.Ny/2) + (0 if (a.Ny % 2) else .5)) * py

orig_xy = np.stack(np.meshgrid(orig_x, orig_y), axis = -1)


omega = ap.omega

bspec = BaseSpec(lMax = a.lMax)
Tmatrix = ap.tmgen(bspec, ap.omega)
particles= [Particle(orig_xy[i], Tmatrix) for i in np.ndindex(orig_xy.shape[:-1])]

sym = FinitePointGroup(point_group_info['D2h'])
ss = ScatteringSystem(particles, sym)

wavenumber = ap.background_epsmu.k(omega).real # Currently, ScatteringSystem does not "remember" frequency nor wavenumber

sinalpha_list = np.linspace(a.kx_lim[0],a.kx_lim[1],a.N)

# Plane wave data
E_cart_list = np.empty((a.N,3), dtype=complex)
E_cart_list[:,:] = np.array((0,1,0))[None,:]
k_cart_list = np.empty((a.N,3), dtype=float)
k_cart_list[:,0] = sinalpha_list
k_cart_list[:,1] = 0
k_cart_list[:,2] = np.sqrt(1-sinalpha_list**2)
k_cart_list *= wavenumber

σ_ext_list_ir = np.empty((a.N, ss.nirreps), dtype=float)
σ_scat_list_ir = np.empty((a.N, ss.nirreps), dtype=float)

outfile_tmp = defaultprefix + ".tmp" if a.output is None else a.output + ".tmp"

for iri in range(ss.nirreps):
    logging.info("processing irrep %d/%d" % (iri, ss.nirreps))
    LU = None # to trigger garbage collection before the next call
    translation_matrix = None
    LU = ss.scatter_solver(wavenumber,iri)
    logging.info("LU solver created")
    translation_matrix = ss.translation_matrix_packed(wavenumber, iri, BesselType.REGULAR) + np.eye(ss.saecv_sizes[iri]) 
    logging.info("auxillary translation matrix created")

    for j in range(a.N):
        # the following two could be calculated only once, but probably not a big deal
        ã = ss.planewave_full(k_cart=k_cart_list[j], E_cart=E_cart_list[j])
        Tã = ss.apply_Tmatrices_full(ã)

        Tãi = ss.pack_vector(Tã, iri)
        ãi = ss.pack_vector(ã, iri)
        fi = LU(Tãi)
        σ_ext_list_ir[j, iri] = -np.vdot(ãi, fi).real/wavenumber**2
        σ_scat_list_ir[j, iri] = np.vdot(fi,np.dot(translation_matrix, fi)).real/wavenumber**2
    if a.save_gradually:
        iriout = outfile_tmp + ".%d" % iri
        np.savez(iriout, iri=iri, meta=vars(a), sinalpha=sinalpha_list, k_cart = k_cart_list, E_cart=E_cart_list,
		 omega=omega, wavenumber=wavenumber, σ_ext_list_ir=σ_ext_list_ir[:,iri], σ_scat_list_ir=σ_scat_list_ir[:,iri])
        logging.info("partial results saved to %s"%iriout)

σ_abs_list_ir = σ_ext_list_ir - σ_scat_list_ir
σ_abs= np.sum(σ_abs_list_ir, axis=-1)
σ_scat= np.sum(σ_scat_list_ir, axis=-1)
σ_ext= np.sum(σ_ext_list_ir, axis=-1)


outfile = defaultprefix + ".npz" if a.output is None else a.output
np.savez(outfile, meta=vars(a), sinalpha=sinalpha_list, k_cart = k_cart_list, E_cart=E_cart_list, σ_ext=σ_ext,σ_abs=σ_abs,σ_scat=σ_scat,
 σ_ext_ir=σ_ext_list_ir,σ_abs_ir=σ_abs_list_ir,σ_scat_ir=σ_scat_list_ir, omega=omega, wavenumber=wavenumber
       )
logging.info("Saved to %s" % outfile)


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

