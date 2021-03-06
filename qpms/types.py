"""
Here shall be defined some types that are needed across the individual 
modules of the qpms package and they do not belong logically to a single
module.
"""

import collections
import enum

class NormalizationT(enum.IntEnum):
    """ Corresponding to the c type qpms_normalization_t from translations.h """
    KRISTENSSON=2 # power-normalised, this should be the default in the future
    TAYLOR=1
    UNDEF=0

class BesselT(enum.IntEnum):
    """ Corresponding to the c type qpms_bessel_t from translations.h """
    BESSEL_REGULAR = 1
    BESSEL_SINGULAR = 2
    HANKEL_PLUS = 3
    HANKEL_MINUS = 4
    UNDEF = 0


'''
The namedtuples below might become classes or other objects in later versions.
'''

TMatrixOp = collections.namedtuple('TMatrixOp',
        ['optype', 'content'])

TMatrixSpec = collections.namedtuple('TMatrixSpec', 
        ['lMax_override', 'tmatrix_path', 'ops'])
TMatrixSpec.__doc__ = """\
Specification of a to-be-created TMatrix object.

lMax_override: int or None. If int and lower than the cutoff degree of the
T-Matrix file located at tmatrix_path, lMax_override shall be used instead.

tmatrix_path: str. Location of the .TMatrix file to be loaded.

ops: sequence of TMatrixOp instances. The operations to be performed on
top of the loaded .TMatrix files.
"""

ParticleSpec = collections.namedtuple('ParticleSpec', ['label', 'position', 'tmatrix_spec'])
ParticleSpec.__doc___ = """\
Specification of an individual scatterer, or a component scatterer 
of a lattice unit cell.

label: immutable (usually a string or None). Unique label of a unit cell
component (and the corresponding sublattice).

position: tuple of floats (or similar, such as numpy array).

tmatrix_spec: TMatrixSpec or TMatrix instance.
"""

LatticeSpec = collections.namedtuple('LatticeSpec', ['basis', 'particle_set_specs'])
LatticeSpec.__doc__ = """\
Specification of a lattice, finite or infinite.

basis: tuple of basic vectors (tuples of floats or similar) (or similar, 
such as 2d numpy array), preferably of reduced basis.

particle_set_specs: sequence of tuples (particle_spec, selection) where
particle_spec is a ParticleSpec and selection is a sequence of 2d integer
coordinates/indices indicating which particles (TODO describe better)
of that type shall be included in a finite lattice. Alternatively, 
selection can be None for an infinite lattice.
"""

