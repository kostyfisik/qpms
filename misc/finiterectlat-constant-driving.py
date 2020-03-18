#!/usr/bin/env python3

import math
from qpms.argproc import ArgParser
figscale=2

ap = ArgParser(['rectlattice2d_finite', 'single_particle', 'single_lMax', 'single_omega'])
ap.add_argument("-k", '--wavevector', nargs=2, type=float, required=True, help='"Bloch" vector, modulating phase of the driving', metavar=('KX', 'KY'), default=(0., 0.))
# ap.add_argument("--kpi", action='store_true', help="Indicates that the k vector is given in natural units instead of SI, i.e. the arguments given by -k shall be automatically multiplied by pi / period (given by -p argument)")
ap.add_argument("-o", "--output", type=str, required=False, help='output path (if not provided, will be generated automatically)')
ap.add_argument("-O", "--plot-out", type=str, required=False, help="path to plot output (optional)")
ap.add_argument("-P", "--plot", action='store_true', help="if -p not given, plot to a default path")
ap.add_argument("-g", "--save-gradually", action='store_true', help="saves the partial result after computing each irrep")

#ap.add_argument("--irrep", type=str, default="none", help="Irrep subspace (irrep index from 0 to 7, irrep label, or 'none' for no irrep decomposition")


a=ap.parse_args()

import logging
logging.basicConfig(format='%(asctime)s %(message)s', level=logging.INFO)

Nx, Ny = a.size
px, py = a.period

particlestr = ("sph" if a.height is None else "cyl") + ("_r%gnm" % (a.radius*1e9))
if a.height is not None: particlestr += "_h%gnm" % (a.height * 1e9)
defaultprefix = "cd_%s_p%gnmx%gnm_%dx%d_m%s_n%g_k_%g_%g_f%geV_L%d" % (
    particlestr, px*1e9, py*1e9, Nx, Ny, str(a.material), a.refractive_index, a.wavevector[0], a.wavevector[1], a.eV, a.lMax,)
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


omega = ap.omega

bspec = BaseSpec(lMax = a.lMax)
medium = EpsMuGenerator(ap.background_epsmu)
particles= [Particle(orig_xy[i], ap.tmgen, bspec) for i in np.ndindex(orig_xy.shape[:-1])]

sym = FinitePointGroup(point_group_info['D2h'])
ss, ssw = ScatteringSystem.create(particles=particles, medium=medium, omega=omega, sym=sym)

wavenumber = ap.background_epsmu.k(omega) # Currently, ScatteringSystem does not "remember" frequency nor wavenumber


outfile_tmp = defaultprefix + ".tmp" if a.output is None else a.output + ".tmp"

nelem = len(bspec)
phases = np.exp(1j*np.dot(ss.positions[:,:2], np.array(a.wavevector)))
driving_full = np.zeros((nelem, ss.fecv_size),dtype=complex)
for y in range(nelem):
    driving_full[y,y::nelem] = phases


scattered_full = np.zeros((nelem, ss.fecv_size),dtype=complex)
scattered_ir = [None for iri in range(ss.nirreps)]


for iri in range(ss.nirreps):
    logging.info("processing irrep %d/%d" % (iri, ss.nirreps))
    LU = None # to trigger garbage collection before the next call
    translation_matrix = None
    LU = ssw.scatter_solver(iri)
    logging.info("LU solver created")
    #translation_matrix = ss.translation_matrix_packed(wavenumber, iri, BesselType.REGULAR) + np.eye(ss.saecv_sizes[iri]) 
    #logging.info("auxillary translation matrix created")

    scattered_ir[iri] = np.empty((nelem, ss.saecv_sizes[iri]), dtype=complex)
    scattered_ir_unpacked = np.empty((nelem, ss.fecv_size), dtype=complex)

    for y in range(nelem):
        ã = driving_full[y]
        Tã  = ssw.apply_Tmatrices_full(ã)
        Tãi = ss.pack_vector(Tã, iri)
        ãi = ss.pack_vector(ã, iri)
        fi = LU(Tãi)
        scattered_ir[iri][y] = fi
        scattered_ir_unpacked[y] = ss.unpack_vector(fi, iri)
        scattered_full[y] += scattered_ir_unpacked[y]
    if a.save_gradually:
        iriout = outfile_tmp + ".%d" % iri
        np.savez(iriout, iri=iri, meta=vars(a), 
		 omega=omega, wavenumber=wavenumber, nelem=nelem, wavevector=np.array(a.wavevector), phases=phases, 
                 positions = ss.positions[:,:2],
                 scattered_ir_packed = scattered_ir[iri], 
                 scattered_ir_full = scattered_ir_unpacked, 
                 )
        logging.info("partial results saved to %s"%iriout)


outfile = defaultprefix + ".npz" if a.output is None else a.output
np.savez(outfile, meta=vars(a), 
		 omega=omega, wavenumber=wavenumber, nelem=nelem, wavevector=np.array(a.wavevector), phases=phases, 
                 positions = ss.positions[:,:2],
                 scattered_ir_packed = scattered_ir, 
                 scattered_full = scattered_full, 
       )
logging.info("Saved to %s" % outfile)


if a.plot or (a.plot_out is not None):
    positions = ss.positions
    xpositions = np.unique(positions[:,0])
    assert(len(xpositions) == Nx)
    ypositions = np.unique(positions[:,1])
    assert(len(ypositions == Ny))
    # particle positions as integer indices
    posmap = np.empty((positions.shape[0],2), dtype=int)
    for i, pos in enumerate(positions):
        posmap[i,0] = np.searchsorted(xpositions, positions[i,0])
        posmap[i,1] = np.searchsorted(ypositions, positions[i,1])

    def fullvec2grid(fullvec):
        arr = np.empty((Nx,Ny,nelem), dtype=complex)
        for pi, offset in enumerate(ss.fullvec_poffsets):
            ix, iy = posmap[pi]
            arr[ix, iy] = fullvec[offset:offset+nelem]
        return arr

    import matplotlib
    matplotlib.use('pdf')
    from matplotlib import pyplot as plt, cm
    t, l, m = bspec.tlm()
    phasecm = cm.twilight

    fig, axes = plt.subplots(nelem, 12, figsize=(figscale*12, figscale*nelem))
    for yp in range(0,3):
        axes[0,4*yp+0].set_title("abs / %s,%d,%+d"%('E' if t[yp]==2 else 'M', l[yp], m[yp],))
        axes[0,4*yp+1].set_title("arg / %s,%d,%+d"%('E' if t[yp]==2 else 'M', l[yp], m[yp],))
        axes[0,4*yp+2].set_title("Fabs / %s,%d,%+d"%('E' if t[yp]==2 else 'M', l[yp], m[yp],))
        axes[0,4*yp+3].set_title("Farg / %s,%d,%+d"%('E' if t[yp]==2 else 'M', l[yp], m[yp],))

    for y in range(nelem):
        axes[y,0].set_ylabel("%s,%d,%+d"%('E' if t[y]==2 else 'M', l[y], m[y],))
        fulvec = scattered_full[y]
        vecgrid = fullvec2grid(fulvec)
        vecgrid_ff = np.fft.fftshift(np.fft.fft2(vecgrid, axes=(0,1)),axes=(0,1))
        lemax = np.amax(abs(vecgrid))
        for yp in range(0,3):
            if(np.amax(abs(vecgrid[...,yp])) > lemax*1e-5):
                axes[y,yp*4].imshow(abs(vecgrid[...,yp]), vmin=0)
                axes[y,yp*4].text(0.5, 0.5, '%g' % np.amax(abs(vecgrid[...,yp])), horizontalalignment='center', verticalalignment='center', transform=axes[y,yp*4].transAxes)
                axes[y,yp*4+1].imshow(np.angle(vecgrid[...,yp]), vmin=-np.pi, vmax=np.pi, cmap=phasecm)
                axes[y,yp*4+2].imshow(abs(vecgrid_ff[...,yp]), vmin=0)
                axes[y,yp*4+3].imshow(np.angle(vecgrid_ff[...,yp]), vmin=-np.pi, vmax=np.pi, cmap=phasecm)
            else:
                for c in range(0,4):
                    axes[y,yp*4+c].tick_params(bottom=False, left=False, labelbottom=False, labelleft=False)
    
    plotfile = defaultprefix + ".pdf" if a.plot_out is None else a.plot_out
    fig.savefig(plotfile)

exit(0)

