#!/usr/bin/env python3

import math
from qpms.argproc import ArgParser


ap = ArgParser(['rectlattice2d_finite', 'single_particle', 'single_lMax', 'single_omega', 'planewave'])
ap.add_argument("-o", "--output", type=str, required=False, help='output path (if not provided, will be generated automatically)')
ap.add_argument("-O", "--plot-out", type=str, required=False, help="path to plot output (optional)")
ap.add_argument("-P", "--plot", action='store_true', help="if -p not given, plot to a default path")
ap.add_argument("-g", "--save-gradually", action='store_true', help="saves the partial result after computing each irrep")


a=ap.parse_args()

import logging
logging.basicConfig(format='%(asctime)s %(message)s', level=logging.INFO)

Nx, Ny = a.size
px, py = a.period


import numpy as np
import qpms
import math
from qpms.qpms_p import cart2sph, sph2cart, sph_loccart2cart, sph_loccart_basis
from qpms.cybspec import BaseSpec
from qpms.cytmatrices import CTMatrix, TMatrixGenerator
from qpms.qpms_c import Particle
from qpms.cymaterials import EpsMu, EpsMuGenerator, LorentzDrudeModel, lorentz_drude
from qpms.cycommon import DebugFlags, dbgmsg_enable
from qpms import FinitePointGroup, ScatteringSystem, BesselType, eV, hbar
from qpms.symmetries import point_group_info
eh = eV/hbar
pi = math.pi

particlestr = ("sph" if a.height is None else "cyl") + ("_r%gnm" % (a.radius*1e9))
if a.height is not None: particlestr += "_h%gnm" % (a.height * 1e9)
defaultprefix = "%s_p%gnmx%gnm_%dx%d_m%s_n%g_φ%gπ_θ(%g_%g)π_ψ%gπ_χ%gπ_f%geV_L%d" % (
    particlestr, px*1e9, py*1e9, Nx, Ny, str(a.material), a.refractive_index, a.phi/pi, np.amin(a.theta)/pi, np.amax(a.theta)/pi, a.psi/pi, a.chi/pi, a.eV, a.lMax, )
logging.info("Default file prefix: %s" % defaultprefix)

dbgmsg_enable(DebugFlags.INTEGRATION)

#Particle positions
orig_x = (np.arange(Nx/2) + (0 if (Nx % 2) else .5)) * px
orig_y = (np.arange(Ny/2) + (0 if (Ny % 2) else .5)) * py

orig_xy = np.stack(np.meshgrid(orig_x, orig_y), axis = -1)


omega = ap.omega

bspec = BaseSpec(lMax = a.lMax)
Tmatrix = ap.tmgen(bspec, ap.omega)
particles= [Particle(orig_xy[i], Tmatrix) for i in np.ndindex(orig_xy.shape[:-1])]

sym = FinitePointGroup(point_group_info['D2h'])
ss, ssw = ScatteringSystem.create(particles, ap.background_emg, omega, sym=sym)

wavenumber = ap.background_epsmu.k(omega).real # Currently, ScatteringSystem does not "remember" frequency nor wavenumber

## Plane wave data
a.theta = np.array(a.theta)
k_sph_list = np.stack((np.broadcast_to(wavenumber, a.theta.shape), a.theta, np.broadcast_to(a.phi, a.theta.shape)), axis=-1)
sψ, cψ = math.sin(a.psi), math.cos(a.psi)
sχ, cχ = math.sin(a.chi), math.cos(a.chi)
E_sph = (0., cψ*cχ + 1j*sψ*sχ, sψ*cχ + 1j*cψ*sχ) 

k_cart_list = sph2cart(k_sph_list)
E_cart_list = sph_loccart2cart(E_sph, k_sph_list)

npoints = a.theta.shape[0]

σ_ext_list_ir = np.empty((npoints, ss.nirreps), dtype=float)
σ_scat_list_ir = np.empty((npoints, ss.nirreps), dtype=float)

outfile_tmp = defaultprefix + ".tmp" if a.output is None else a.output + ".tmp"

for iri in range(ss.nirreps):
    logging.info("processing irrep %d/%d" % (iri, ss.nirreps))
    LU = None # to trigger garbage collection before the next call
    translation_matrix = None
    LU = ssw.scatter_solver(iri)
    logging.info("LU solver created")
    translation_matrix = ss.translation_matrix_packed(wavenumber, iri, BesselType.REGULAR) + np.eye(ss.saecv_sizes[iri]) 
    logging.info("auxillary translation matrix created")

    for j in range(npoints):
        # the following two could be calculated only once, but probably not a big deal
        ã = ss.planewave_full(k_cart=k_cart_list[j], E_cart=E_cart_list[j])
        Tã = ssw.apply_Tmatrices_full(ã)

        Tãi = ss.pack_vector(Tã, iri)
        ãi = ss.pack_vector(ã, iri)
        fi = LU(Tãi)
        σ_ext_list_ir[j, iri] = -np.vdot(ãi, fi).real/wavenumber**2
        σ_scat_list_ir[j, iri] = np.vdot(fi,np.dot(translation_matrix, fi)).real/wavenumber**2
    if a.save_gradually:
        iriout = outfile_tmp + ".%d" % iri
        np.savez(iriout, iri=iri, meta=vars(a), k_sph=k_sph_list, k_cart = k_cart_list, E_cart=E_cart_list, E_sph=np.array(E_sph),
		 omega=omega, wavenumber=wavenumber, σ_ext_list_ir=σ_ext_list_ir[:,iri], σ_scat_list_ir=σ_scat_list_ir[:,iri])
        logging.info("partial results saved to %s"%iriout)

σ_abs_list_ir = σ_ext_list_ir - σ_scat_list_ir
σ_abs= np.sum(σ_abs_list_ir, axis=-1)
σ_scat= np.sum(σ_scat_list_ir, axis=-1)
σ_ext= np.sum(σ_ext_list_ir, axis=-1)


outfile = defaultprefix + ".npz" if a.output is None else a.output
np.savez(outfile, meta=vars(a), k_sph=k_sph_list, k_cart = k_cart_list, E_cart=E_cart_list, E_sph=np.array(E_sph),
        σ_ext=σ_ext,σ_abs=σ_abs,σ_scat=σ_scat,
        σ_ext_ir=σ_ext_list_ir,σ_abs_ir=σ_abs_list_ir,σ_scat_ir=σ_scat_list_ir, omega=omega, wavenumber=wavenumber
       )
logging.info("Saved to %s" % outfile)


if a.plot or (a.plot_out is not None):
    import matplotlib
    matplotlib.use('pdf')
    from matplotlib import pyplot as plt

    fig = plt.figure()
    ax = fig.add_subplot(111)
    sintheta = np.sin(a.theta)
    ax.plot(sintheta, σ_ext*1e12,label='$\sigma_\mathrm{ext}$')
    ax.plot(sintheta, σ_scat*1e12, label='$\sigma_\mathrm{scat}$')
    ax.plot(sintheta, σ_abs*1e12, label='$\sigma_\mathrm{abs}$')
    ax.legend()
    ax.set_xlabel('$\sin\\theta$')
    ax.set_ylabel('$\sigma/\mathrm{\mu m^2}$')
    
    plotfile = defaultprefix + ".pdf" if a.plot_out is None else a.plot_out
    fig.savefig(plotfile)

exit(0)

