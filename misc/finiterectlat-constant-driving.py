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
ap.add_argument("-S", "--symmetry-adapted", default=None, help="Use a symmetry-adapted basis of a given point group instead of individual spherical harmonics")
ap.add_argument("-d", "--ccd-distance", type=float, default=math.nan, help='Far-field "CCD" distance from the sample')
ap.add_argument("-D", "--ccd-size", type=float, default=math.nan, help='Far-field "CCD" width and heighth')
ap.add_argument("-R", "--ccd-resolution", type=int, default=101, help='Far-field "CCD" resolution')

#ap.add_argument("--irrep", type=str, default="none", help="Irrep subspace (irrep index from 0 to 7, irrep label, or 'none' for no irrep decomposition")


a=ap.parse_args()

import logging
logging.basicConfig(format='%(asctime)s %(message)s', level=logging.INFO)

Nx, Ny = a.size
px, py = a.period

particlestr = ("sph" if a.height is None else "cyl") + ("_r%gnm" % (a.radius*1e9))
if a.height is not None: particlestr += "_h%gnm" % (a.height * 1e9)
defaultprefix = "cd_%s_p%gnmx%gnm_%dx%d_m%s_n%g_k_%g_%g_f%geV_L%d_micro-%s" % (
    particlestr, px*1e9, py*1e9, Nx, Ny, str(a.material), a.refractive_index, a.wavevector[0], a.wavevector[1], a.eV, a.lMax, "SO3" if a.symmetry_adapted is None else a.symmetry_adapted)
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

def float_nicestr(x, tol=1e-5):
    x = float(x)
    if .5**2 - abs(x) < tol:
        return(("-" if x < 0 else '+') + "2^{-2}")
    else: 
        return "%+.3g" % x

def cplx_nicestr(x):
    x = complex(x)
    if x == 0:
        return '0'
    ret = ""
    if x.real:
        ret = ret + float_nicestr(x.real)
    if x.imag:
        ret = ret + float_nicestr(x.imag) + 'i'
    if x.real and x.imag:
        return '(' + ret + ')'
    else:
        return ret

def cleanarray(a, atol=1e-10, copy=True):
    a = np.array(a, copy=copy)
    sieve = abs(a.real) < atol
    a[sieve] = 1j * a[sieve].imag
    sieve = abs(a.imag) < atol
    a[sieve] = a[sieve].real
    return a

def nicerot(a, atol=1e-10, copy=True): #gives array a "nice" phase
    a = np.array(a, copy=copy)
    i = np.argmax(abs(a))
    a = a / a[i] * abs(a[i])
    return a

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
if a.symmetry_adapted is not None:
    ss1, ssw1 = ScatteringSystem.create(particles=[Particle((0,0,0), ap.tmgen, bspec)], medium=medium, omega=omega, 
                    sym=FinitePointGroup(point_group_info[a.symmetry_adapted]))
    fvcs1 = np.empty((nelem, nelem), dtype=complex)
    y = 0
    iris1 = []
    for iri1 in range(ss1.nirreps):
        for j in range(ss1.saecv_sizes[iri1]):
            pvc1 = np.zeros((ss1.saecv_sizes[iri1],), dtype=complex)
            pvc1[j] = 1
            fvcs1[y] = ss1.unpack_vector(pvc1, iri1)
            fvcs1[y] = cleanarray(nicerot(fvcs1[y], copy=False), copy=False)
            driving_full[y] = (phases[:, None] * fvcs1[y][None,:]).flatten()
            y += 1
            iris1.append(iri1)
    iris1 = np.array(iris1)
else:
    for y in range(nelem):
        driving_full[y,y::nelem] = phases


scattered_full = np.zeros((nelem, ss.fecv_size),dtype=complex)
scattered_ir = [None for iri in range(ss.nirreps)]

ir_contained = np.ones((nelem, ss.nirreps), dtype=bool)

for iri in range(ss.nirreps):
    logging.info("processing irrep %d/%d" % (iri, ss.nirreps))
    LU = None # to trigger garbage collection before the next call
    translation_matrix = None
    LU = ssw.scatter_solver(iri)
    logging.info("LU solver created")
    #translation_matrix = ss.translation_matrix_packed(wavenumber, iri, BesselType.REGULAR) + np.eye(ss.saecv_sizes[iri]) 
    #logging.info("auxillary translation matrix created")

    scattered_ir[iri] = np.zeros((nelem, ss.saecv_sizes[iri]), dtype=complex)
    scattered_ir_unpacked = np.zeros((nelem, ss.fecv_size), dtype=complex)

    for y in range(nelem):
        ã = driving_full[y]
        ãi = cleanarray(ss.pack_vector(ã, iri), copy=False)
        if np.all(ãi == 0):
            ir_contained[y, iri] = False
            continue
        Tã  = ssw.apply_Tmatrices_full(ã)
        Tãi = ss.pack_vector(Tã, iri)
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

t, l, m = bspec.tlm()

if not math.isnan(a.ccd_distance):
    logging.info("Computing the far fields")
    ccd_size = (50 * a.ccd_distance / (max(Nx*px, Ny*py) * ssw.wavenumber.real)) if math.isnan(a.ccd_size) else a.ccd_size
    ccd_x = np.linspace(-ccd_size/2, ccd_size/2, a.ccd_resolution)
    ccd_y = np.linspace(-ccd_size/2, ccd_size/2, a.ccd_resolution)
    ccd_grid = np.meshgrid(ccd_x, ccd_y, (a.ccd_distance,), indexing='ij')
    ccd_points = np.stack(ccd_grid, axis=-1).squeeze(axis=-2)
    ccd_fields = np.empty((nelem,) + ccd_points.shape, dtype=complex)
    for y in range(nelem):
        ccd_fields[y] = ssw.scattered_E(scattered_full[y], ccd_points, btyp=BesselType.HANKEL_PLUS)
    logging.info("Far fields done")

outfile = defaultprefix + ".npz" if a.output is None else a.output
np.savez(outfile, meta=vars(a), 
		 omega=omega, wavenumber=wavenumber, nelem=nelem, wavevector=np.array(a.wavevector), phases=phases, 
                 positions = ss.positions[:,:2],
                 scattered_ir_packed = scattered_ir, 
                 scattered_full = scattered_full,
                 ir_contained = ir_contained,
                 t=t, l=l, m=m,
                 iris1 = iris1 if (a.symmetry_adapted is not None) else None,
                 irnames1 = ss1.irrep_names if (a.symmetry_adapted is not None) else None,
                 fvcs1 = fvcs1 if (a.symmetry_adapted is not None) else None,
                 #ccd_size = ccd_size if not math.isnan(a.ccd_distance) else None,
                 ccd_points = ccd_points if not math.isnan(a.ccd_distance) else None,
                 ccd_fields = ccd_fields if not math.isnan(a.ccd_distance) else None,
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
    pmcm = cm.bwr
    abscm = cm.plasma 

    fig, axes = plt.subplots(nelem, 12 if math.isnan(a.ccd_distance) else 15, figsize=(figscale*(12 if math.isnan(a.ccd_distance) else 15), figscale*nelem))
    for yp in range(0,3): # TODO xy-dipoles instead?
        axes[0,4*yp+0].set_title("abs / %s,%d,%+d"%('E' if t[yp]==2 else 'M', l[yp], m[yp],))
        axes[0,4*yp+1].set_title("arg / %s,%d,%+d"%('E' if t[yp]==2 else 'M', l[yp], m[yp],))
        axes[0,4*yp+2].set_title("Fabs / %s,%d,%+d"%('E' if t[yp]==2 else 'M', l[yp], m[yp],))
        axes[0,4*yp+3].set_title("Farg / %s,%d,%+d"%('E' if t[yp]==2 else 'M', l[yp], m[yp],))
    if not math.isnan(a.ccd_distance):
        #axes[0,12].set_title("$E_{xy}$ @ $z = %g; \phi$" % a.ccd_distance)
        #axes[0,13].set_title("$E_{xy}$ @ $z = %g; \phi + \pi/2$" % a.ccd_distance)
        axes[0,12].set_title("$|E_{x}|^2$ @ $z = %g\,\mathrm{m}$" % a.ccd_distance)
        axes[0,13].set_title("$|E_{y}|^2$ @ $z = %g\,\mathrm{m}$" % a.ccd_distance)
        axes[0,14].set_title("$|E_{z}|^2$ @ $z = %g\,\mathrm{m}$" % a.ccd_distance)

    for y in range(nelem):
        fulvec = scattered_full[y]
        if a.symmetry_adapted is not None:
            driving_nonzero_y = [j for j in range(nelem) if abs(fvcs1[y,j]) > 1e-5]
            driving_descr = ss1.irrep_names[iris1[y]]+'\n'+', '.join(('$'+cplx_nicestr(fvcs1[y,j])+'$' +
"(%s,%d,%+d)" % (("E" if t[j] == 2 else "M"), l[j], m[j]) for j in
driving_nonzero_y)) # TODO shorten the complex number precision
        else:
            driving_descr = "%s,%d,%+d"%('E' if t[y]==2 else 'M', l[y], m[y],)
        axes[y,0].set_ylabel(driving_descr)
        vecgrid = fullvec2grid(fulvec)
        vecgrid_ff = np.fft.fftshift(np.fft.fft2(vecgrid, axes=(0,1)),axes=(0,1))
        lemax = np.amax(abs(vecgrid))
        for yp in range(0,3):
            if(np.amax(abs(vecgrid[...,yp])) > lemax*1e-5):
                axes[y,yp*4].imshow(abs(vecgrid[...,yp]), vmin=0, interpolation='none')
                axes[y,yp*4].text(0.5, 0.5, '%g' % np.amax(abs(vecgrid[...,yp])), horizontalalignment='center', verticalalignment='center', transform=axes[y,yp*4].transAxes)
                axes[y,yp*4+1].imshow(np.angle(vecgrid[...,yp]), vmin=-np.pi, vmax=np.pi, cmap=phasecm, interpolation='none')
                axes[y,yp*4+2].imshow(abs(vecgrid_ff[...,yp]), vmin=0, interpolation='none')
                axes[y,yp*4+3].imshow(np.angle(vecgrid_ff[...,yp]), vmin=-np.pi, vmax=np.pi, cmap=phasecm, interpolation='none')
            else:
                for c in range(0,4):
                    axes[y,yp*4+c].tick_params(bottom=False, left=False, labelbottom=False, labelleft=False)
        if not math.isnan(a.ccd_distance):
            fxye=(-ccd_size/2, ccd_size/2, -ccd_size/2, ccd_size/2)
            e2vmax = np.amax(np.linalg.norm(ccd_fields[y], axis=-1)**2)
            xint = abs(ccd_fields[y,...,0])**2
            yint = abs(ccd_fields[y,...,1])**2
            axes[y, 12].imshow(xint, origin="lower",vmax=e2vmax, extent=fxye, cmap=abscm, interpolation='none')
            axes[y, 13].imshow(yint, origin="lower",vmax=e2vmax, extent=fxye, cmap=abscm, interpolation='none')
            #axes[y, 12].imshow(np.sum(abs(ccd_fields[y,...,:2].real)**2, axis=-1), origin="lower",vmax=e2vmax, extent=fxye, cmap=abscm)
            #axes[y, 12].quiver(ccd_points[...,1], ccd_points[...,0], ccd_fields[y,...,1].real, ccd_fields[y,...,0].real, color='w')
            #axes[y, 13].imshow(np.sum(abs(ccd_fields[y,...,:2].imag)**2, axis=-1) ,origin="lower",vmax=e2vmax, extent=fxye, cmap=abscm)
            #axes[y, 13].quiver(ccd_points[...,1], ccd_points[...,0],ccd_fields[y,...,1].imag, ccd_fields[y,...,0].imag, color='w')
            zint = abs(ccd_fields[y,...,2])**2
            axes[y, 14].imshow(zint, origin='lower', extent=fxye, cmap=abscm, interpolation='none')
            axes[y, 14].text(0.5, 0.5, '%g' % (np.amax(zint)/e2vmax), 
                    horizontalalignment='center', verticalalignment='center', transform=axes[y,14].transAxes)
    plotfile = defaultprefix + ".pdf" if a.plot_out is None else a.plot_out
    fig.savefig(plotfile)

exit(0)

