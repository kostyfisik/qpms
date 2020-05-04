#!/usr/bin/env python3

import math
from qpms.argproc import ArgParser, make_dict_action, sslice
figscale=3

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
ap.add_argument("--xslice", default={None:None}, nargs=2, 
        action=make_dict_action(argtype=sslice, postaction='append', first_is_key=True), 
        )
ap.add_argument("--yslice", default={None:None}, nargs=2, 
        action=make_dict_action(argtype=sslice, postaction='append', first_is_key=True), 
        )


#ap.add_argument("--irrep", type=str, default="none", help="Irrep subspace (irrep index from 0 to 7, irrep label, or 'none' for no irrep decomposition")


a=ap.parse_args()

import logging
logging.basicConfig(format='%(asctime)s %(message)s', level=logging.INFO)

Nx, Ny = a.size
px, py = a.period

particlestr = ("sph" if a.height is None else "cyl") + ("_r%gnm" % (a.radius*1e9))
if a.height is not None: particlestr += "_h%gnm" % (a.height * 1e9)
defaultprefix = "cd_%s_p%gnmx%gnm_%dx%d_m%s_n%s_k_%g_%g_f%geV_L%d_micro-%s" % (
    particlestr, px*1e9, py*1e9, Nx, Ny, str(a.material), str(a.background), a.wavevector[0], a.wavevector[1], a.eV, a.lMax, "SO3" if a.symmetry_adapted is None else a.symmetry_adapted)
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

# Check slice ranges and generate all corresponding combinations
slicepairs = []
slicelabels = set(a.xslice.keys())  | set(a.yslice.keys())
for label in slicelabels:
    rowslices = a.xslice.get(label, None)
    colslices = a.yslice.get(label, None)
    # TODO check validity of the slices.
    if rowslices is None:
        rowslices = [slice(None, None, None)]
    if colslices is None:
        colslices = [slice(None, None, None)]
    for rs in rowslices:
        for cs in colslices:
            slicepairs.append((rs, cs))

def realdipfieldlabels(yp):
    if yp == 0: return 'x'
    if yp == 1: return 'y'
    if yp == 2: return 'z'
    raise ValueError
def realdipfields(vecgrid, yp):
    if yp == 1:
        return vecgrid[...,0] + vecgrid[...,2]
    if yp == 0:
        return -1j*(vecgrid[...,0] - vecgrid[...,2])
    if yp == 2:
        return vecgrid[...,1]
    raise ValueError

def float_nicestr(x, tol=1e-5):
    x = float(x)
    if .5**2 - abs(x) < tol:
        return(("-" if x < 0 else '+') + "2^{-2}")
    else: 
        return "%+.3g" % x

def cplx_nicestr(x, tol=1e-5):
    x = complex(x)
    if x == 0:
        return '0'
    ret = ""
    if x.real:
        ret = ret + float_nicestr(x.real, tol)
    if x.imag:
        ret = ret + float_nicestr(x.imag, tol) + 'i'
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

# Mapping between ss particles and grid positions
positions = ss.positions
xpositions = np.unique(positions[:,0])
assert(len(xpositions) == Nx)
ypositions = np.unique(positions[:,1])
assert(len(ypositions == Ny))
# particle positions as integer indices
posmap = np.empty((positions.shape[0],2), dtype=int)
invposmap = np.empty((Nx, Ny), dtype=int)
for i, pos in enumerate(positions):
    posmap[i,0] = np.searchsorted(xpositions, positions[i,0])
    posmap[i,1] = np.searchsorted(ypositions, positions[i,1])
    invposmap[posmap[i,0], posmap[i, 1]] = i

def fullvec2grid(fullvec, swapxy=False):
    arr = np.empty((Nx,Ny,nelem), dtype=complex)
    for pi, offset in enumerate(ss.fullvec_poffsets):
        ix, iy = posmap[pi]
        arr[ix, iy] = fullvec[offset:offset+nelem]
    return np.swapaxes(arr, 0, 1) if swapxy else arr


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


# Apply the driving on the specified slices only
nsp = len(slicepairs)
driving_full_sliced = np.zeros((nsp,) + driving_full.shape, dtype=complex)
p1range = np.arange(nelem)
for spi in range(nsp):
    xs, ys = slicepairs[spi]
    driven_pi = invposmap[xs, ys].flatten()
    driven_y = ((driven_pi * nelem)[:,None] + p1range[None,:]).flatten()
    driving_full_sliced[spi][:, driven_y] = driving_full[:, driven_y]

scattered_full = np.zeros((nsp, nelem, ss.fecv_size),dtype=complex)
scattered_ir = [None for iri in range(ss.nirreps)]

ir_contained = np.ones((nsp, nelem, ss.nirreps), dtype=bool)

for iri in range(ss.nirreps):
    logging.info("processing irrep %d/%d" % (iri, ss.nirreps))
    LU = None # to trigger garbage collection before the next call
    translation_matrix = None
    LU = ssw.scatter_solver(iri)
    logging.info("LU solver created")
    #translation_matrix = ss.translation_matrix_packed(wavenumber, iri, BesselType.REGULAR) + np.eye(ss.saecv_sizes[iri]) 
    #logging.info("auxillary translation matrix created")

    scattered_ir[iri] = np.zeros((nsp, nelem, ss.saecv_sizes[iri]), dtype=complex)
    scattered_ir_unpacked = np.zeros((nsp, nelem, ss.fecv_size), dtype=complex)

    for spi in range(nsp):
        for y in range(nelem):
            ã = driving_full_sliced[spi,y]
            ãi = cleanarray(ss.pack_vector(ã, iri), copy=False)
            if np.all(ãi == 0):
                ir_contained[spi, y, iri] = False
                continue
            Tã  = ssw.apply_Tmatrices_full(ã)
            Tãi = ss.pack_vector(Tã, iri)
            fi = LU(Tãi)
            scattered_ir[iri][spi, y] = fi
            scattered_ir_unpacked[spi, y] = ss.unpack_vector(fi, iri)
            scattered_full[spi, y] += scattered_ir_unpacked[spi, y]
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
    ccd_points = np.swapaxes(np.stack(ccd_grid, axis=-1).squeeze(axis=-2), 0,1) # First axis is y, second is x, because of imshow...
    ccd_fields = np.empty((nsp, nelem,) + ccd_points.shape, dtype=complex)
    for spi in range(nsp):
        for y in range(nelem):
            ccd_fields[spi, y] = ssw.scattered_E(scattered_full[spi, y], ccd_points, btyp=BesselType.HANKEL_PLUS)
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


    import matplotlib
    matplotlib.use('pdf')
    from matplotlib import pyplot as plt, cm
    from matplotlib.backends.backend_pdf import PdfPages
    t, l, m = bspec.tlm()
    phasecm = cm.twilight
    pmcm = cm.bwr
    abscm = cm.plasma 
    
    plotfile = defaultprefix + ".pdf" if a.plot_out is None else a.plot_out
    pp = PdfPages(plotfile)

    for spi in range(nsp):
        fig, axes = plt.subplots(nelem, 12 if math.isnan(a.ccd_distance) else 16, figsize=(figscale*(12 if math.isnan(a.ccd_distance) else 16), figscale*nelem))
        for yp in range(0,3): # TODO xy-dipoles instead?
            axes[0,4*yp+0].set_title("abs / (E,1,%s)" % realdipfieldlabels(yp))
            axes[0,4*yp+1].set_title("arg / (E,1,%s)" % realdipfieldlabels(yp))
            axes[0,4*yp+2].set_title("Fabs / (E,1,%s)" % realdipfieldlabels(yp))
            axes[0,4*yp+3].set_title("Farg / (E,1,%s)" % realdipfieldlabels(yp))
        if not math.isnan(a.ccd_distance):
            #axes[0,12].set_title("$E_{xy}$ @ $z = %g; \phi$" % a.ccd_distance)
            #axes[0,13].set_title("$E_{xy}$ @ $z = %g; \phi + \pi/2$" % a.ccd_distance)
            axes[0,12].set_title("$|E_{x}|^2$ @ $z = %g\,\mathrm{m}$" % a.ccd_distance)
            axes[0,13].set_title("$|E_{y}|^2$ @ $z = %g\,\mathrm{m}$" % a.ccd_distance)
            axes[0,14].set_title("$|E_x + E_y|^2$ @ $z = %g\,\mathrm{m}$" % a.ccd_distance)
            axes[0,15].set_title("$|E_{z}|^2$ @ $z = %g\,\mathrm{m}$" % a.ccd_distance)
            for gg in range(12,16):
                axes[-1,gg].set_xlabel("$x/\mathrm{m}$")
            

        for y in range(nelem):
            fulvec = scattered_full[spi,y]
            if a.symmetry_adapted is not None:
                driving_nonzero_y = [j for j in range(nelem) if abs(fvcs1[y,j]) > 1e-5]
                driving_descr = ss1.irrep_names[iris1[y]]+'\n'+', '.join(('$'+cplx_nicestr(fvcs1[y,j])+'$' +
    "(%s,%d,%+d)" % (("E" if t[j] == 2 else "M"), l[j], m[j]) for j in
    driving_nonzero_y)) # TODO shorten the complex number precision
            else:
                driving_descr = "%s,%d,%+d"%('E' if t[y]==2 else 'M', l[y], m[y],)
            axes[y,0].set_ylabel(driving_descr)
            axes[y,-1].yaxis.set_label_position("right")
            axes[y,-1].set_ylabel("$y/\mathrm{m}$\n"+driving_descr)
            vecgrid = fullvec2grid(fulvec, swapxy=True)
            vecgrid_ff = np.fft.fftshift(np.fft.fft2(vecgrid, axes=(0,1)),axes=(0,1))
            lemax = np.amax(abs(vecgrid))
            for yp in range(0,3):
                if(np.amax(abs(realdipfields(vecgrid,yp))) > lemax*1e-5):
                    axes[y,yp*4].imshow(abs(realdipfields(vecgrid,yp)), vmin=0, interpolation='none')
                    axes[y,yp*4].text(0.5, 0.5, '%g' % np.amax(abs(realdipfields(vecgrid,yp))), horizontalalignment='center', verticalalignment='center', transform=axes[y,yp*4].transAxes)
                    axes[y,yp*4+1].imshow(np.angle(realdipfields(vecgrid,yp)), vmin=-np.pi, vmax=np.pi, cmap=phasecm, interpolation='none')
                    axes[y,yp*4+2].imshow(abs(realdipfields(vecgrid_ff,yp)), vmin=0, interpolation='none')
                    axes[y,yp*4+3].imshow(np.angle(realdipfields(vecgrid_ff,yp)), vmin=-np.pi, vmax=np.pi, cmap=phasecm, interpolation='none')
                else:
                    for c in range(0,4):
                        axes[y,yp*4+c].tick_params(bottom=False, left=False, labelbottom=False, labelleft=False)
            if not math.isnan(a.ccd_distance):
                fxye=(-ccd_size/2, ccd_size/2, -ccd_size/2, ccd_size/2)
                e2vmax = np.amax(np.linalg.norm(ccd_fields[spi,y], axis=-1)**2)
                xint = abs(ccd_fields[spi,y,...,0])**2
                yint = abs(ccd_fields[spi,y,...,1])**2
                xyint = abs(ccd_fields[spi,y,...,0] + ccd_fields[spi,y,...,1])**2
                zint = abs(ccd_fields[spi,y,...,2])**2
                xintmax = np.amax(xint)
                yintmax = np.amax(yint)
                zintmax = np.amax(zint)
                xyintmax = np.amax(xyint)
                axes[y, 12].imshow(xint, origin="lower", extent=fxye, cmap=abscm, interpolation='none')
                axes[y, 13].imshow(yint, origin="lower", extent=fxye, cmap=abscm, interpolation='none')
                axes[y, 14].imshow(xyint, origin="lower", extent=fxye, cmap=abscm, interpolation='none')
                axes[y, 15].imshow(zint, origin='lower', extent=fxye, cmap=abscm, interpolation='none')
                axes[y, 12].text(0.5, 0.5, '%g\n%g' % (xintmax,xintmax/e2vmax), 
                        horizontalalignment='center', verticalalignment='center', transform=axes[y,12].transAxes)
                axes[y, 13].text(0.5, 0.5, '%g\n%g' % (yintmax,yintmax/e2vmax), 
                        horizontalalignment='center', verticalalignment='center', transform=axes[y,13].transAxes)
                axes[y, 14].text(0.5, 0.5, '%g\n%g' % (xyintmax,xyintmax/e2vmax), 
                        horizontalalignment='center', verticalalignment='center', transform=axes[y,14].transAxes)
                axes[y, 15].text(0.5, 0.5, '%g\n%g' % (zintmax,zintmax/e2vmax), 
                        horizontalalignment='center', verticalalignment='center', transform=axes[y,15].transAxes)
                for gg in range(12,16):
                    axes[y,gg].yaxis.tick_right()
                for gg in range(12,15):
                    axes[y,gg].yaxis.set_major_formatter(plt.NullFormatter())
        fig.text(0, 0, str(slicepairs[spi]), horizontalalignment='left', verticalalignment='bottom')
        pp.savefig()
    pp.close()

exit(0)

