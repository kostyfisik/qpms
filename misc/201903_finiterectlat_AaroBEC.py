#!/usr/bin/env python
# coding: utf-8
from qpms import Particle, CTMatrix, BaseSpec, FinitePointGroup, ScatteringSystem, TMatrixInterpolator, eV, hbar, c
from qpms.symmetries import point_group_info
import numpy as np
import os
import sys
nm = 1e-9

sym = FinitePointGroup(point_group_info['D2h'])
bspec = BaseSpec(lMax = 2)
tmfile = '/m/phys/project/qd/Marek/tmatrix-experiments/Cylinder/AaroBEC/cylinder_50nm_lMax4_cleaned.TMatrix'
#outputdatadir = '/home/necadam1/wrkdir/AaroBECfinite_new'
outputdatadir = '/u/46/necadam1/unix/project/AaroBECfinite_new'
os.makedirs(outputdatadir, exist_ok = True)
interp = TMatrixInterpolator(tmfile, bspec, symmetrise = sym, atol = 1e-8)
# There is only one t-matrix in the system for each frequency. We initialize the matrix with the lowest frequency data.
# Later, we can replace it using the tmatrix[...] = interp(freq) and s.update_tmatrices NOT YET; TODO

omega = float(sys.argv[3]) * eV/hbar
sv_threshold = float(sys.argv[4])

# Now place the particles and set background index.
px = 571*nm; py = 621*nm
n = 1.51
Nx = int(sys.argv[1])
Ny = int(sys.argv[2])

orig_x = (np.arange(Nx/2) + (0 if (Nx % 2) else .5)) * px
orig_y = (np.arange(Ny/2) + (0 if (Ny % 2) else .5)) * py

orig_xy = np.stack(np.meshgrid(orig_x, orig_y), axis = -1)

tmatrix = interp(omega)
particles = [Particle(orig_xy[i], tmatrix) for i in np.ndindex(orig_xy.shape[:-1])]


ss = ScatteringSystem(particles, sym)


k = n * omega / c


for iri in range(ss.nirreps):
    mm_iri = ss.modeproblem_matrix_packed(k, iri)
    U, S, Vh = np.linalg.svd(mm_iri)
    print(iri, ss.irrep_names[iri], S[-1])
    starti = max(0,len(S) - np.searchsorted(S[::-1], sv_threshold, side='left')-1)
    np.savez(os.path.join(outputdatadir, 'Nx%d_Ny%d_%geV_ir%d.npz'%(Nx, Ny, omega/eV*hbar, iri)),
        S=S[starti:], omega=omega, Vh = Vh[starti:], iri=iri, Nx = Nx, Ny= Ny )
    # Don't forget to conjugate Vh before transforming it to the full vector!

