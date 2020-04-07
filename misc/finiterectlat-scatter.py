#!/usr/bin/env python3

import math
pi = math.pi
from qpms.argproc import ArgParser


ap = ArgParser(['rectlattice2d_finite', 'single_particle', 'single_lMax', 'omega_seq_real_ng', 'planewave'])
ap.add_argument("-o", "--output", type=str, required=False, help='output path (if not provided, will be generated automatically)')
ap.add_argument("-O", "--plot-out", type=str, required=False, help="path to plot output (optional)")
ap.add_argument("-P", "--plot", action='store_true', help="if -p not given, plot to a default path")
ap.add_argument("-g", "--save-gradually", action='store_true', help="saves the partial result after computing each irrep")


a=ap.parse_args()

import logging
logging.basicConfig(format='%(asctime)s %(message)s', level=logging.INFO)

import numpy as np
import qpms
from qpms.qpms_p import cart2sph, sph2cart, sph_loccart2cart, sph_loccart_basis
from qpms.cybspec import BaseSpec
from qpms.cytmatrices import CTMatrix, TMatrixGenerator
from qpms.qpms_c import Particle
from qpms.cymaterials import EpsMu, EpsMuGenerator, LorentzDrudeModel, lorentz_drude
from qpms.cycommon import DebugFlags, dbgmsg_enable
from qpms import FinitePointGroup, ScatteringSystem, BesselType, eV, hbar
from qpms.symmetries import point_group_info
eh = eV/hbar

dbgmsg_enable(DebugFlags.INTEGRATION)

Nx, Ny = a.size
px, py = a.period

particlestr = ("sph" if a.height is None else "cyl") + ("_r%gnm" % (a.radius*1e9))
if a.height is not None: particlestr += "_h%gnm" % (a.height * 1e9)
defaultprefix = "%s_p%gnmx%gnm_%dx%d_m%s_bg%s_φ%gπ_θ(%g_%g)π_ψ%gπ_χ%gπ_f%s_L%d" % (
    particlestr, px*1e9, py*1e9, Nx, Ny, str(a.material), str(a.background), a.phi/pi, np.amin(a.theta)/pi, np.amax(a.theta)/pi, a.psi/pi, a.chi/pi, ap.omega_descr, a.lMax, )
logging.info("Default file prefix: %s" % defaultprefix)


#Particle positions
orig_x = (np.arange(Nx/2) + (0 if (Nx % 2) else .5)) * px
orig_y = (np.arange(Ny/2) + (0 if (Ny % 2) else .5)) * py

orig_xy = np.stack(np.meshgrid(orig_x, orig_y), axis = -1)

bspec = BaseSpec(lMax = a.lMax)
particles= [Particle(orig_xy[i], ap.tmgen, bspec=bspec) for i in np.ndindex(orig_xy.shape[:-1])]

sym = FinitePointGroup(point_group_info['D2h'])
ss, ssw = ScatteringSystem.create(particles, ap.background_emg, ap.allomegas[0], sym=sym)


## Plane wave data
a.theta = np.atleast_1d(np.array(a.theta))
dir_sph_list = np.stack((np.broadcast_to(1, a.theta.shape), a.theta, np.broadcast_to(a.phi, a.theta.shape)), axis=-1)
sψ, cψ = math.sin(a.psi), math.cos(a.psi)
sχ, cχ = math.sin(a.chi), math.cos(a.chi)
E_sph = (0., cψ*cχ + 1j*sψ*sχ, sψ*cχ + 1j*cψ*sχ) 

dir_cart_list = sph2cart(dir_sph_list)
E_cart_list = sph_loccart2cart(E_sph, dir_sph_list)

nfreq = len(ap.allomegas)
ndir = a.theta.shape[0]

k_cart_arr = np.empty((nfreq, ndir, 3), dtype=float)
wavenumbers = np.empty((nfreq,), dtype=float)

σ_ext_arr_ir = np.empty((nfreq, ndir, ss.nirreps), dtype=float)
σ_scat_arr_ir = np.empty((nfreq, ndir, ss.nirreps), dtype=float)

outfile_tmp = defaultprefix + ".tmp" if a.output is None else a.output + ".tmp"

for i, omega in enumerate(ap.allomegas):
    logging.info("Processing frequency %g eV" % (omega / eV,))
    if i != 0:
        ssw = ss(omega)
    if ssw.wavenumber.imag != 0:
        warnings.warn("The background medium wavenumber has non-zero imaginary part. Don't expect emaningful results for cross sections.")
    wavenumber = ssw.wavenumber.real
    wavenumbers[i] = wavenumber

    k_sph_list = np.array(dir_sph_list, copy=True)
    k_sph_list[:,0] = wavenumber

    k_cart_arr[i] = sph2cart(k_sph_list)

    for iri in range(ss.nirreps):
        logging.info("processing irrep %d/%d" % (iri, ss.nirreps))
        LU = None # to trigger garbage collection before the next call
        translation_matrix = None
        LU = ssw.scatter_solver(iri)
        logging.info("LU solver created")
        translation_matrix = ssw.translation_matrix_packed(iri, BesselType.REGULAR) + np.eye(ss.saecv_sizes[iri]) 
        logging.info("auxillary translation matrix created")

        for j in range(ndir):
            k_cart = k_cart_arr[i,j]
            # the following two could be calculated only once, but probably not a big deal
            ã = ss.planewave_full(k_cart=k_cart_arr[i,j], E_cart=E_cart_list[j])
            Tã = ssw.apply_Tmatrices_full(ã)

            Tãi = ss.pack_vector(Tã, iri)
            ãi = ss.pack_vector(ã, iri)
            fi = LU(Tãi)
            σ_ext_arr_ir[i, j, iri] = -np.vdot(ãi, fi).real/wavenumber**2
            σ_scat_arr_ir[i, j, iri] = np.vdot(fi,np.dot(translation_matrix, fi)).real/wavenumber**2
        if a.save_gradually:
            iriout = outfile_tmp + ".%d.%d" % (i, iri)
            np.savez(iriout, omegai=i, iri=iri, meta=vars(a), omega=omega, k_sph=k_sph_list, k_cart = k_cart_arr, E_cart=E_cart_list, E_sph=np.array(E_sph),
                     wavenumber=wavenumber, σ_ext_list_ir=σ_ext_arr_ir[i,:,iri], σ_scat_list_ir=σ_scat_list_ir[i,:,iri])
            logging.info("partial results saved to %s"%iriout)

σ_abs_arr_ir = σ_ext_arr_ir - σ_scat_arr_ir
σ_abs_arr = np.sum(σ_abs_arr_ir, axis=-1)
σ_scat_arr = np.sum(σ_scat_arr_ir, axis=-1)
σ_ext_arr = np.sum(σ_ext_arr_ir, axis=-1)


outfile = defaultprefix + ".npz" if a.output is None else a.output
np.savez(outfile, meta=vars(a), k_sph=k_sph_list, k_cart = k_cart_arr, E_cart=E_cart_list, E_sph=np.array(E_sph),
        σ_ext=σ_ext_arr,σ_abs=σ_abs_arr,σ_scat=σ_scat_arr,
        σ_ext_ir=σ_ext_arr_ir,σ_abs_ir=σ_abs_arr_ir,σ_scat_ir=σ_scat_arr_ir, omega=ap.allomegas, wavenumbers=wavenumbers
       )
logging.info("Saved to %s" % outfile)

if a.plot or (a.plot_out is not None):
    import matplotlib
    from matplotlib.backends.backend_pdf import PdfPages
    matplotlib.use('pdf')
    from matplotlib import pyplot as plt
    from scipy.interpolate import griddata

    plotfile = defaultprefix + ".pdf" if a.plot_out is None else a.plot_out
    with PdfPages(plotfile) as pdf:
        ipm = 'nearest'
        sintheta = np.sin(a.theta)
        if False: #len(ap.omega_ranges) != 0:
            # angle plot ---------------------------------
            fig = plt.figure(figsize=(210/25.4, 297/25.4))
            vmax = max(np.amax(σ_ext_arr), np.amax(σ_scat_arr), np.amax(σ_abs_arr))
            vmin = min(np.amin(σ_ext_arr), np.amin(σ_scat_arr), np.amin(σ_abs_arr))

            ax = fig.add_subplot(311)
            ax.pcolormesh(a.theta, ap.allomegas/eh, σ_ext_arr, vmin=vmin, vmax=vmax)
            ax.set_xlabel('$\\theta$')
            ax.set_ylabel('$\\hbar\\omega / \\mathrm{eV}$')
            ax.set_title('$\\sigma_\\mathrm{ext}$')

            ax = fig.add_subplot(312)
            ax.pcolormesh(a.theta, ap.allomegas/eh, σ_scat_arr, vmin=vmin, vmax=vmax)
            ax.set_xlabel('$\\theta$')
            ax.set_ylabel('$\\hbar\\omega / \\mathrm{eV}$')
            ax.set_title('$\\sigma_\\mathrm{scat}$')

            ax = fig.add_subplot(313)
            im = ax.pcolormesh(a.theta, ap.allomegas/eh, σ_abs_arr, vmin=vmin, vmax=vmax)
            ax.set_xlabel('$\\theta$')
            ax.set_ylabel('$\\hbar\\omega / \\mathrm{eV}$')
            ax.set_title('$\\sigma_\\mathrm{abs}$')


            fig.subplots_adjust(right=0.8)
            fig.colorbar(im, cax = fig.add_axes([0.85, 0.15, 0.05, 0.7]))

            pdf.savefig(fig)
            plt.close(fig)

        if len(ap.omega_ranges) != 0:
            # "k-space" plot -----------------------------
            domega = np.amin(np.diff(ap.allomegas))
            dsintheta = np.amin(abs(np.diff(sintheta)))
            dk = dsintheta * wavenumbers[0]

            # target image grid
            grid_y, grid_x = np.mgrid[ap.allomegas[0] : ap.allomegas[-1] : domega, np.amin(sintheta) * wavenumbers[-1] : np.amax(sintheta) * wavenumbers[-1] : dk]
            imextent = (np.amin(sintheta) * wavenumbers[-1] / 1e6, np.amax(sintheta) * wavenumbers[-1] / 1e6, ap.allomegas[0] / eh, ap.allomegas[-1] / eh)

            # source coordinates for griddata
            ktheta = sintheta[None, :] * wavenumbers[:, None]
            omegapoints = np.broadcast_to(ap.allomegas[:, None], ktheta.shape)
            points = np.stack( (ktheta.flatten(), omegapoints.flatten()), axis = -1)

            fig = plt.figure(figsize=(210/25.4, 297/25.4))
            vmax = np.amax(σ_ext_arr)

            ax = fig.add_subplot(311)
            grid_z = griddata(points, σ_ext_arr.flatten(), (grid_x, grid_y), method = ipm)
            ax.imshow(grid_z, extent = imextent, origin = 'lower', vmin = 0, vmax = vmax, aspect = 'auto', interpolation='none')
            ax.set_xlabel('$k_\\theta / \\mathrm{\\mu m^{-1}}$')
            ax.set_ylabel('$\\hbar\\omega / \\mathrm{eV}$')
            ax.set_title('$\\sigma_\\mathrm{ext}$')

            ax = fig.add_subplot(312)
            grid_z = griddata(points, σ_scat_arr.flatten(), (grid_x, grid_y), method = ipm)
            ax.imshow(grid_z, extent = imextent, origin = 'lower', vmin = 0, vmax = vmax, aspect = 'auto', interpolation='none')
            ax.set_xlabel('$k_\\theta / \\mathrm{\\mu m^{-1}}$')
            ax.set_ylabel('$\\hbar\\omega / \\mathrm{eV}$')
            ax.set_title('$\\sigma_\\mathrm{scat}$')

            ax = fig.add_subplot(313)
            grid_z = griddata(points, σ_abs_arr.flatten(), (grid_x, grid_y), method = ipm)
            im = ax.imshow(grid_z, extent = imextent, origin = 'lower', vmin = 0, vmax = vmax, aspect = 'auto', interpolation='none')
            ax.set_xlabel('$k_\\theta / \\mathrm{\\mu m^{-1}}$')
            ax.set_ylabel('$\\hbar\\omega / \\mathrm{eV}$')
            ax.set_title('$\\sigma_\\mathrm{abs}$')

            fig.subplots_adjust(right=0.8)
            fig.colorbar(im, cax = fig.add_axes([0.85, 0.15, 0.05, 0.7]))

            pdf.savefig(fig)
            plt.close(fig)

        for omega in ap.omega_singles:
            i = np.searchsorted(ap.allomegas, omega)

            fig = plt.figure()
            fig.suptitle("%g eV" % (omega / eh))
            ax = fig.add_subplot(111)
            sintheta = np.sin(a.theta)
            ax.plot(sintheta, σ_ext_arr[i]*1e12,label='$\sigma_\mathrm{ext}$')
            ax.plot(sintheta, σ_scat_arr[i]*1e12, label='$\sigma_\mathrm{scat}$')
            ax.plot(sintheta, σ_abs_arr[i]*1e12, label='$\sigma_\mathrm{abs}$')
            ax.legend()
            ax.set_xlabel('$\sin\\theta$')
            ax.set_ylabel('$\sigma/\mathrm{\mu m^2}$')

            pdf.savefig(fig)
            plt.close(fig)


exit(0)

