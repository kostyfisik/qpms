from qpms import Particle, CTMatrix, BaseSpec, FinitePointGroup, ScatteringSystem, TMatrixInterpolator, eV, hbar, c
from qpms.symmetries import point_group_info
import numpy as np
import os
import sys
nm = 1e-9

sym = FinitePointGroup(point_group_info['D2h'])
bspec = BaseSpec(lMax = 2)
#tmfile = '/m/phys/project/qd/Marek/tmatrix-experiments/Cylinder/AaroBEC/cylinder_50nm_lMax4_cleaned.TMatrix'
tmfile = '/home/mmn/repo/tmatrix-experiments/Cylinder/AaroBEC/cylinder_50nm_lMax4_cleaned.TMatrix'
#outputdatadir = '/home/necadam1/wrkdir/AaroBECfinite_new'
outputdatadir = '/tmp/u/46/necadam1/unix/project/AaroBECfinite_new'
os.makedirs(outputdatadir, exist_ok = True)
interp = TMatrixInterpolator(tmfile, bspec, symmetrise = sym, atol = 1e-8)
# There is only one t-matrix in the system for each frequency. We initialize the matrix with the lowest frequency data.
# Later, we can replace it using the tmatrix[...] = interp(freq) and s.update_tmatrices NOT YET; TODO

omega = 1.475 * eV/hbar
sv_threshold = 0.5

# Now place the particles and set background index.
px = 571*nm; py = 621*nm
n = 1.51
Nx = 5
Ny = 7

orig_x = (np.arange(Nx/2) + (0 if (Nx % 2) else .5)) * px
orig_y = (np.arange(Ny/2) + (0 if (Ny % 2) else .5)) * py

orig_xy = np.stack(np.meshgrid(orig_x, orig_y), axis = -1)

tmatrix = interp(omega)
particles = [Particle(orig_xy[i], tmatrix) for i in np.ndindex(orig_xy.shape[:-1])]


ss = ScatteringSystem(particles, sym)
k = n * omega / c

for iri in range(ss.nirreps):
    mm_iri_orig = ss.modeproblem_matrix_packed(k, iri, version = None)
    mm_iri_alt = ss.modeproblem_matrix_packed(k, iri, version='R')
    mm_iri_paral = ss.modeproblem_matrix_packed(k, iri, version='pR')
    print(np.amax(abs(mm_iri_orig-mm_iri_alt)), np.amax(abs(mm_iri_orig-mm_iri_paral)))

