'''
Object oriented approach for the classical multiple scattering problem.
'''
import numpy as np
import time
import scipy
import sys
from qpms_c import * # TODO be explicit about what is imported
from .qpms_p import nelem2lMax # TODO be explicit about what is imported

class Scattering(object):
    '''

    This is the most general class for a system of scatterers
    in a non-lossy homogeneous background
    to be solved with the multiple_scattering method. The scatterers,
    as long as they comply with the disjoint circumscribed sphere
    hypothesis, can each have any position in the 3D space and
    any T-matrix.

    Note that this object describes the scattering problem only for
    a single given frequency, as the T-matrices and wavelenght
    otherwise differ and all the computationally demanding
    parts have to be done for each frequency. However,
    the object can be recycled for many incident field shapes
    at the given frequency.

    Attributes should be perhaps later redefined to be read-only
    (or make descriptors for them).

    Args:
        positions: (N,3)-shaped real array
        TMatrices: (N,2,nelem,2,nelem)-shaped array
        k_0 (float): Wave number for the space between scatterers.

    Attributes:
        positions:
        TMatrices:
        k_0 (float): Wave number for the space between scatterers.
        lMax (int): Absolute maximum l for all scatterers. Depending on implementation,
            lMax can be smaller for some individual scatterers in certain subclasses.
        prepared (bool): Keeps information whether the interaction matrix has
            already been built and factorized.
        

    '''
    def __init__(self, positions, TMatrices, k_0, lMax = None, verbose=False):
        self.positions = positions
        self.TMatrices = TMatrices
        self.k_0 = k_0
        self.lMax = lMax ? lMax : nelem2lMax(TMatrices.shape[-1])
        self.prepared = False

    def prepare(self, keep_interaction_matrix = False, verbose=False):
        if not self.prepared:
            if not self.interaction_matrix:
                self.build_interaction_matrix(verbose=verbose)
            self.lupiv = scipy.linalg_lu_factor(interaction_matrix)
            if not keep_interaction_matrix:
                self.interaction_matrix = None
            self.prepared = True

    def build_interaction_matrix(verbose = False):
        pass

    def scatter(pq_0_c, verbose = False):
        self.prepare(verbose=verbose)


        pass




