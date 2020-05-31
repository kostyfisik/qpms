TODO list before 1.0 release
============================

- Tests!
- Docs!
- Cross section calculations. (Done in some Python scripts.)
- Field calculations. (Partly done, needs more testing.)
  * Also test periodic vs. nonperiodic consistence (big finite lattice + absorbing medium vs. infinite lattice + absorbing medium).
- Complex frequencies, n's, k's. (Mostly done.)
- Transforming point (meta)generators.
- Check whether moble's quaternions and my 
  quaternions give the same results in tmatrices.py
- Ewald summations of all types of lattices (dimensionality-wise).
- Split lattices.h into separate point generator and lattice vector manipulation parts.
  * Maybe move something from the .h to .c file.
- Check exact normalisation convention of scuff-tmatrix output.
- Check whether the Condon-Shortley phase affects the form of Wigner matrices.
- The xflip, yflip and possible i-factor problem.
- General 3D point group symmetries.
  * Instead the current hard-coded limited set.
  * The generation, finding subgroups etc. should be "easy" with
    quaternions and stuff, as the  set is quite limited, 
    see [Wikipedia](https://en.wikipedia.org/wiki/Point_groups_in_three_dimensions).
  * Not sure about the representations, though.
  * As a description of a T-matrix / particle metadata.
- Nice CLI for all general enough utilities.
- Remove legacy code.
- Split qpms_c.pyx.
- Reduce compiler warnings.
- Python exceptions instead of hard crashes in the C library where possible.
- Scatsystem init sometimes fail due to rounding errors and hardcoded absolute tolerance 
  in the qpms_tmatrix_isclose() call.
- Prefix all identifiers. Maybe think about a different prefix than qpms?
- Consistent indentation and style overall.
- Rewrite the parallelized translation matrix, mode problem matrix generators
  in a way that reuses as much code as possible without copypasting

Nice but less important features
--------------------------------

- Static, thread-safe caches of constant coefficients + API without the current "calculators".

