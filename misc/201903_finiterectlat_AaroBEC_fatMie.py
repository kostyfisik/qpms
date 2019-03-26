#!/usr/bin/env python3
# coding: utf-8
from qpms import Particle, CTMatrix, BaseSpec, FinitePointGroup, ScatteringSystem, TMatrixInterpolator, eV, hbar, c, MaterialInterpolator, scatsystem_set_nthreads
from qpms.symmetries import point_group_info
import numpy as np
import os
import sys
from pathlib import Path

if 'SLURM_CPUS_PER_TASK' in os.environ:
    scatsystem_set_nthreads(os.environ['SLURM_CPUS_PER_TASK'])

nm = 1e-9

rewrite_output = '--rewrite-output' in sys.argv

radiusfactor = float(sys.argv[5])

cyr_part_height = 50*nm
cyr_part_radius = 50*nm
cyr_part_volume = cyr_part_height * np.pi * cyr_part_radius**2
eqv_sph_radius = (3/4/np.pi*cyr_part_volume)**(1/3) * radiusfactor

sym = FinitePointGroup(point_group_info['D2h'])
bspec = BaseSpec(lMax = 2)
#tmfile = '/m/phys/project/qd/Marek/tmatrix-experiments/Cylinder/AaroBEC/cylinder_50nm_lMax4_cleaned.TMatrix'
materialfile = '/home/necadam1/wrkdir/repo/refractiveindex.info-database/database/data/main/Au/Johnson.yml'

#outputdatadir = '/home/necadam1/wrkdir/AaroBECfinite_new'
#outputdatadir = '/u/46/necadam1/unix/project/AaroBECfinite_sph'
outputdatadir = '/home/necadam1/wrkdir/AaroBECfinite_fatsph'
os.makedirs(outputdatadir, exist_ok = True)
mi = MaterialInterpolator(materialfile)
#interp = TMatrixInterpolator(tmfile, bspec, symmetrise = sym, atol = 1e-8)
# There is only one t-matrix in the system for each frequency. We initialize the matrix with the lowest frequency data.
# Later, we can replace it using the tmatrix[...] = interp(freq) and s.update_tmatrices NOT YET; TODO

omega = float(sys.argv[3]) * eV/hbar
sv_threshold = float(sys.argv[4])


# Now place the particles and set background index.
px = 571*nm; py = 621*nm
n = 1.52
Nx = int(sys.argv[1])
Ny = int(sys.argv[2])

orig_x = (np.arange(Nx/2) + (0 if (Nx % 2) else .5)) * px
orig_y = (np.arange(Ny/2) + (0 if (Ny % 2) else .5)) * py

orig_xy = np.stack(np.meshgrid(orig_x, orig_y), axis = -1)

#tmatrix = interp(omega)
tmatrix = CTMatrix.spherical_perm(bspec, eqv_sph_radius, omega, mi(omega), n**2)
particles = [Particle(orig_xy[i], tmatrix) for i in np.ndindex(orig_xy.shape[:-1])]


ss = ScatteringSystem(particles, sym)


k = n * omega / c


for iri in range(ss.nirreps):
    destpath = os.path.join(outputdatadir, 'Nx%d_Ny%d_%geV_ir%d_r%gnm.npz'%(Nx, Ny, omega/eV*hbar, iri, eqv_sph_radius/nm))
    touchpath = os.path.join(outputdatadir, 'Nx%d_Ny%d_%geV_ir%d_r%gnm.done'%(Nx, Ny, omega/eV*hbar, iri, eqv_sph_radius/nm))
    if (os.path.isfile(destpath) or os.path.isfile(touchpath)) and not rewrite_output:
      print(destpath, 'already exists, skipping')
      continue
    mm_iri = ss.modeproblem_matrix_packed(k, iri)
    U, S, Vh = np.linalg.svd(mm_iri)
    del U
    print(iri, ss.irrep_names[iri], S[-1])
    starti = max(0,len(S) - np.searchsorted(S[::-1], sv_threshold, side='left')-1)
    np.savez(destpath,
        S=S[starti:], omega=omega, Vh = Vh[starti:], iri=iri, Nx = Nx, Ny= Ny )
    del S
    del Vh
    Path(touchpath).touch()
    # Don't forget to conjugate Vh before transforming it to the full vector!

