from qpms_cdefs cimport *

from .cybspec cimport BaseSpec
from .cybspec import VSWFNorm, default_bspec
from .cycommon import VSWFType, BesselType

def vswf_single(kr_csph, qpms_vswf_type_t t, qpms_l_t l, qpms_m_t m, 
        qpms_bessel_t btyp = QPMS_BESSEL_REGULAR, qpms_normalisation_t norm = default_bspec.norm):
    cdef csph_t kr
    kr.r, kr.theta, kr.phi = kr_csph
    if t == QPMS_VSWF_ELECTRIC:
        return qpms_vswf_single_el_csph(m, l, kr, btyp, norm)
    elif t == QPMS_VSWF_MAGNETIC:
        return qpms_vswf_single_mg_csph(m, l, kr, btyp, norm)
    elif t == QPMS_VSWF_LONGITUDINAL:
        raise NotImplementedError("Longitudinal single waves not yet implemented, sorry.")
    else:
        raise ValueError("Invalid wave type specified")


