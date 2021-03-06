/*! \file qpms_types.h
 * \brief Common qpms types.
 */
#ifndef QPMS_TYPES_H
#define QPMS_TYPES_H
#include <complex.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#ifndef M_PI_2
#define M_PI_2 (1.570796326794896619231321691639751442098584699687552910487L)
#endif
#ifndef M_PI
#define M_PI (3.14159265358979323846264338327950288419716939937510582097494L)
#endif
#ifndef M_SQRT2
#define M_SQRT2 (1.41421356237309504880168872420969807856967187537694807317668L)
#endif
#ifndef M_SQRTPI
#define M_SQRTPI (1.77245385090551602729816748334114518279754945612238712821381L)
#endif

// integer index types
typedef int qpms_lm_t;
/// Type for spherical harmonic degree \a l.
typedef int qpms_l_t; /* can't be unsigned because of the behaviour under the - operator;
			 also -1 needed as an invalid value for scalar waves. */

/// Type for spherical harmonic order \a m.
typedef qpms_lm_t qpms_m_t;

/// Type for the (\a l, \a m) multiindex of transversal (\a M or \a N -type) VSWFs.
/** This corresponds to the typical memory layout for various coefficient etc.
 *  Corresponds to the l-primary, m-secondary ordering, i.e.
 *  \f[ y = 0: l = 1, m = -1, \f]
 *  \f[ y = 1: l = 1, m =  0, \f]
 *  \f[ y = 2: l = 1, m = +1, \f]
 *  \f[ y = 3: l = 2, m = -2, \f]
 *  ...
 *
 *  See also indexing.h.
 */
typedef size_t qpms_y_t;

/// Type for the (\a l, \a m) multiindex of spherical harmonics, including (0, 0).
/** This differs from qpms_y_t by being shifted by one and including
 *  the \a l = 0 option. Suitable also for scalar and longitudinal waves.
 *  Corresponds to the \a l -primary, \a m -secondary ordering, i.e.
 *  \f[ y = 0: l = 0, m = 0,  \f]
 *  \f[ y = 1: l = 1, m = -1, \f]
 *  \f[ y = 2: l = 1, m = 0,  \f]
 *  \f[ y = 3: l = 1, m = +1, \f]
 *  \f[ y = 4: l = 2, m = -2, \f]
 *  ...
 *
 *  See also indexing.h.
 */
typedef size_t qpms_y_sc_t;

/// Codes of the VSWF types (electric/N, magnetic/M, longitudinal/L).
typedef enum {
	QPMS_VSWF_ELECTRIC = 2, ///< "Electric" (\a N -type) transversal wave.
	QPMS_VSWF_MAGNETIC = 1, ///< "Magnetic" (\a M -type) transversal wave.
	QPMS_VSWF_LONGITUDINAL = 0 ///< Longitudinal (\a L -type) wave (not relevant for radiation).
} qpms_vswf_type_t;

// FIXME DOC The references to the functions do not work in doxygen:
/// Exhaustive index type for VSWF basis functions.
/** Carries information about the wave being of \a M/N/L (magnetic, electric,
 *  or longitudinal) type, as well as the wave's degree and order (\a l, \a m).
 *
 *  The formula is 4 * (qpms_y_sc_t) y_sc + (qmps_vswf_type_t) type_code,
 *  but don't rely on this and use the functions
 *  \ref qpms_tmn2uvswfi() and \ref qpms_uvswfi2tmn()
 *  from qpms_types.h instead
 *  as the formula might change in future versions.
 *
 *  See also indexing.h.
 */
typedef unsigned long long qpms_uvswfi_t; 

/// Error codes / return values for certain numerical functions.
/** Those with values between -2 and 32 are subset of the GSL error codes. */
typedef enum {
	QPMS_SUCCESS = 0, ///< Success.
	QPMS_ERROR = -1, ///< Unspecified error.
	QPMS_FAILURE = QPMS_ERROR,
	QPMS_ENOMEM = 8, ///< Out of memory.
	QPMS_ESING = 21, ///< Apparent singularity detected.
	QPMS_NAN_ENCOUNTERED = 1024 ///< NaN value encountered in data processing.
} qpms_errno_t;



/// Vector spherical wavefuction normalisation and phase convention codes.
/**
 *  Throughout the literature, various conventions for VSWF bases are used.
 *  These bit flags are used by the functions declared in normalisation.h
 *  that return the appropriate convention-dependent factors.
 *
 *  See @ref vswf_conventions for comparison of the various conventions used.
 */
typedef enum {
	QPMS_NORMALISATION_UNDEF = 0, ///< Convention undefined. This should not happen.
       	/// Flag indicating that qpms_normalisition_factor_* should actually return values inverse to the default.
	QPMS_NORMALISATION_INVERSE = 1, 
	/** Flag indicating inversion of the asimuthal phase for complex spherical harmonics (i.e. \f$ e^{-im\phi} \f$
	 * instead of \f$ e^{im\phi} \f$.
	 */
	QPMS_NORMALISATION_REVERSE_AZIMUTHAL_PHASE = 2,
	/// Flag indicating use of the real spherical harmonics.
	/** If QPMS_NORMALISATION_REVERSE_AZIMUTHAL_PHASE is unset, negative \a m
	 * correspond to sine in the asimuthal factor; if set, undefined behaviour.
	 */
	QPMS_NORMALISATION_SPHARM_REAL = 4,
	/// Flag indicating usage of Condon-Shortley phase.
	/** If set, the Ferrers functions and everything derived from them 
	 * (spherical harmonics, VSWFs) will include a \f$ (-1)^m \f$ factor.
	 *
	 * On implementation level, this means that the relevant `gsl_sf_legendre_*_e()`
	 * functions will be called with argument `csphase = -1.` instead of `+1.`.
	 */
	QPMS_NORMALISATION_CSPHASE = 8,
	QPMS_NORMALISATION_M_I = 16, ///< Include an additional \a i -factor into the magnetic waves.
	QPMS_NORMALISATION_M_MINUS = 32, ///< Include an additional \f$-1\f$ -factor into the magnetic waves.
	QPMS_NORMALISATION_N_I = 64, ///< Include an additional \a i -factor into the electric waves.
	QPMS_NORMALISATION_N_MINUS = 128, ///< Include an additional \f$-1\f$ -factor into the magnetic waves.
	QPMS_NORMALISATION_L_I = 256, ///< Include an additional \a i -factor into the longitudinal waves.
	QPMS_NORMALISATION_L_MINUS = 512, ///< Include an additional \f$-1\f$ -factor into the longitudinal waves.
	QPMS_NORMALISATION_NORM_BITSTART = 65536, 
	/// The VSWFs shall be power-normalised. This is the "default".
	/**
	 * Power normalisation is used e.g. in \cite kristensson_spherical_2014 (complex spherical
	 * harmonics with Condon-Shortley phase) or \cite kristensson_scattering_2016 (real
	 * spherical harmonics). This is also the reference for all the other normalisation conventions,
	 * meaning that qpms_normalisation_factor_M() and qpms_normalisation_factor_N() shall
	 * always return `1. + 0.*I` if `norm == QPMS_NORMALISATION_NORM_POWER`.
	 */
	QPMS_NORMALISATION_NORM_POWER = QPMS_NORMALISATION_NORM_BITSTART * 1,
	/// The VSWFs shall be normalised as in \cite taylor_optical_2011 .
	/** This includes a \f$ \sqrt{l(l+1)} \f$ factor compared to the power normalisation. */
	QPMS_NORMALISATION_NORM_SPHARM = QPMS_NORMALISATION_NORM_BITSTART * 3,
	/// The VSWFs shall be created using spherical harmonics without any normalisation. Do not use.
	/** This includes a \f[ 
	 * 	\sqrt{l(l+1)} \left(\frac{(2l+1)}{4\pi}\frac{(l-m)!}{(l+m)!}\right)^{-\frac{1}{2}}
	 * \f] factor compared to the power normalisation.
	 *
	 * Note that this has no sense whatsoever for real spherical harmonics.
	 * Again, do not use this.
	 */
	QPMS_NORMALISATION_NORM_NONE = QPMS_NORMALISATION_NORM_BITSTART * 2,
	QPMS_NORMALISATION_NORM_BITS = QPMS_NORMALISATION_NORM_POWER 
		| QPMS_NORMALISATION_NORM_NONE | QPMS_NORMALISATION_NORM_SPHARM,

	/// VSWF convention used in \cite kristensson_scattering_2016
	QPMS_NORMALISATION_CONVENTION_KRISTENSSON_REAL = QPMS_NORMALISATION_NORM_POWER
		| QPMS_NORMALISATION_SPHARM_REAL,
	/// VSWF convention used in \cite kristensson_spherical_2014
	QPMS_NORMALISATION_CONVENTION_KRISTENSSON = QPMS_NORMALISATION_NORM_POWER
		| QPMS_NORMALISATION_CSPHASE,
	/// VSWF convention used in SCUFF-EM \cite reid_electromagnetism_2016
	QPMS_NORMALISATION_CONVENTION_SCUFF = QPMS_NORMALISATION_NORM_POWER
		| QPMS_NORMALISATION_CSPHASE | QPMS_NORMALISATION_M_I
		| QPMS_NORMALISATION_N_MINUS,
	/// Default VSWF convention. We might encourage the compiler to expect this one.
	QPMS_NORMALISATION_DEFAULT = QPMS_NORMALISATION_CONVENTION_KRISTENSSON
} qpms_normalisation_t;

/// Determine whether the convention includes Condon-Shortley phase (-1) or not (+1).
static inline int qpms_normalisation_t_csphase(qpms_normalisation_t norm) {
	return (norm & QPMS_NORMALISATION_CSPHASE)? -1 : 1;
}

/// Bessel function kinds.
typedef enum {
	QPMS_BESSEL_REGULAR = 1, ///< regular (spherical) Bessel function \a j (Bessel function of the first kind)
	QPMS_BESSEL_SINGULAR = 2, ///< singular (spherical)  Bessel function \a y (Bessel function of the second kind)
	QPMS_HANKEL_PLUS = 3, ///< (spherical) Hankel function \f$ h_1 = j + iy \f$
	QPMS_HANKEL_MINUS = 4, ///< (spherical) Hankel function \f$ h_2 = j - iy \f$
	QPMS_BESSEL_UNDEF = 0 ///< invalid / unspecified kind
} qpms_bessel_t;

// coordinate system types
/// 3D cartesian coordinates. See also vectors.h.
typedef struct cart3_t {
	double x, y, z;
} cart3_t;

/// 3D complex (actually 6D) coordinates. See also vectors.h.
typedef struct ccart3_t {
	complex double x, y, z;
} ccart3_t;

/// 3D complex vector pair (represents the E, H fields).
typedef struct ccart3_pair {
	ccart3_t E, H;
} ccart3_pair;

/// 2D cartesian coordinates. See also vectors.h.
/** See also vectors.h */
typedef struct cart2_t {
	double x, y;
} cart2_t;

/// Spherical coordinates. See also vectors.h.
typedef struct sph_t {
	double r, theta, phi;
} sph_t;

/// Spherical coordinates with complex radial component. See also vectors.h.
typedef struct csph_t { // Do I really need this???
	complex double r;
	double	theta, phi;
} csph_t;

/// 3D complex vector components in local spherical basis. See also vectors.h.
typedef struct csphvec_t {
	complex double rc, thetac, phic; 
} csphvec_t;

/// 2D polar coordinates. See also vectors.h.
typedef struct pol_t {
	double r, phi;
} pol_t;

/// Union type capable to contain various 1D, 2D and 3D coordinates.
/** Usually combined with qpms_coord_system_t. */
typedef union anycoord_point_t {
	double z; ///< 1D cartesian coordinate.
	cart3_t cart3;
	cart2_t cart2;
	sph_t sph;
	pol_t pol;
} anycoord_point_t;

/// Enum codes for common coordinate systems.
typedef enum { 
	// IF EVER CHANGING THE CONSTANT VALUES HERE, 
	// CHECK THAT THEY DO NOT CLASH WITH THOSE IN PGenPointFlags!
	QPMS_COORDS_CART1 = 64, ///< 1D cartesian (= double).
	QPMS_COORDS_POL = 128, ///< 2D polar.
	QPMS_COORDS_SPH = 256, ///< 3D spherical.
	QPMS_COORDS_CART2 = 512, ///< 2D cartesian.
	QPMS_COORDS_CART3 = 1024, ///< 3D cartesian.
	/// Convenience bitmask (not a valid flag!).
	QPMS_COORDS_BITRANGE = QPMS_COORDS_CART1
		| QPMS_COORDS_POL
		| QPMS_COORDS_SPH
		| QPMS_COORDS_CART2
		| QPMS_COORDS_CART3,
} qpms_coord_system_t;

/// Quaternion type.
/**
 * Internaly represented as a pair of complex numbers,
 * \f$ Q_a = Q_1 + iQ_z, Q_b = Q_y + i Q_x\f$.
 *
 * See quaternions.h for "methods".
 */
typedef struct qpms_quat_t {
        complex double a, b;
} qpms_quat_t;

/// Quaternion type as four doubles.
/** See quaternions.h for "methods".
 */
typedef struct qpms_quat4d_t {
        double c1, ci, cj, ck;
} qpms_quat4d_t;

/// 3D improper rotations represented as a quaternion and a sign of the determinant.
/** See quaternions.h for "methods".
 */
typedef struct qpms_irot3_t {
        qpms_quat_t rot; ///< Quaternion representing the rotation part.
        short det; ///< Determinant of the transformation (valid values are 1 (rotation) or -1 (improper rotation)
} qpms_irot3_t;

/// Specifies a finite set of VSWFs.
/**
 * When for example not all the M and N -type waves up to a degree lMax
 * need to be computed, this will specify the subset.
 *
 * A typical use case would be when only even/odd fields wrt. z-plane
 * mirror symmetry are considered.
 * 
 * See vswf.h for "methods".
 */
typedef struct qpms_vswf_set_spec_t {
        size_t n; ///< Actual number of VSWF indices included in ilist.
        qpms_uvswfi_t *ilist; ///< List of wave indices.
	/// Maximum degree of the waves specified in ilist overall.
	/** `max(lMax_M, lMax_N, lMax_L)` */
        qpms_l_t lMax; 
	/// Maximum degree of the magnetic (M-type) waves. 
	/** Set to 0 if no magnetic waves present. */
        qpms_l_t lMax_M, 
	/// Maximum degree of the electric (N-type) waves. 
	/** Set to 0 if no electric waves present. */
                 lMax_N,
	/// Maximum degree of the longitudinal (L-type) waves. 
	/** Set to -1 if no longitudinal waves present. */
                 lMax_L;
        size_t capacity; ///< Allocated capacity of ilist.
        qpms_normalisation_t norm; ///< Normalisation convention. To be set manually if needed.
} qpms_vswf_set_spec_t;

/// T-matrix index used in qpms_scatsys_t and related structures.
typedef int32_t qpms_ss_tmi_t;

/// T-matrix generator index used in qpms_scatsys_t and related structures.
typedef int32_t qpms_ss_tmgi_t;

/// Particle index used in qpms_scatsys_t and related structures.
typedef int32_t qpms_ss_pi_t;

// These types are mainly used in groups.h:
/// Finite group member index. See also groups.h.
typedef int qpms_gmi_t;

/// Irreducible representation index. See also groups.h.
typedef int qpms_iri_t;
/// Constant labeling that no irrep decomposition is done (full system solved instead).
#define QPMS_NO_IRREP ((qpms_iri_t) -1)

/// Permutation representation of a group element.
/** For now, it's just a string of the form "(0,1)(3,4,5)"
 */
typedef const char * qpms_permutation_t;

/// A T-matrix.
/** In the future, I might rather use a more abstract approach in which T-matrix
 *  is a mapping (function) of the field expansion coefficients.
 *  So the interface might change.
 *  For now, let me stick to the square dense matrix representation.
 */
typedef struct qpms_tmatrix_t {
        /** \brief VSWF basis specification, NOT owned by qpms_tmatrix_t by default.
         *
         *  Usually not checked for meaningfulness by the functions (methods),
         *  so the caller should take care that \a spec->ilist does not
         *  contain any duplicities and that for each wave with order \a m
         *  there is also one with order \a −m.
         */
        const qpms_vswf_set_spec_t *spec;
        complex double *m; ///< Matrix elements in row-major order.
        bool owns_m; ///< Information wheter m shall be deallocated with qpms_tmatrix_free()
} qpms_tmatrix_t;


/// Classification of a 3D point group.
typedef enum {

	// Axial groups
	QPMS_PGS_CN, ///< Rotational symmetry \f$ \mathrm{C_{n}} \f$.
	QPMS_PGS_S2N, ///< Rotoreflectional symmetry \f$ \mathrm{S_{2n}} \f$.
	QPMS_PGS_CNH, ///< Rotational symmetry with horizontal reflection \f$ \mathrm{C_{nh}} \f$.
	QPMS_PGS_CNV, ///< Pyramidal symmetry \f$ \mathrm{C_{nv}} \f$.
	QPMS_PGS_DN, ///< Dihedral symmetry \f$ \mathrm{D_{n}} \f$.
	QPMS_PGS_DND, ///< Antiprismatic symmetry \f$ \mathrm{D_{nd}} \f$.
	QPMS_PGS_DNH, ///< Prismatic symmetry \f$ \mathrm{D_{nh}} \f$.

	// Remaining polyhedral groups
	QPMS_PGS_T, ///< Chiral tetrahedral symmetry \f$ \mathrm{T} \f$.
	QPMS_PGS_TD, ///< Full tetrahedral symmetry \f$ \mathrm{T_d} \f$.
	QPMS_PGS_TH, ///< Pyritohedral symmetry \f$ \mathrm{T_h} \f$.
	QPMS_PGS_O,  ///< Chiral octahedral symmetry \f$ \mathrm{O_h} \f$.
	QPMS_PGS_OH,  ///< Full octahedral symmetry \f$ \mathrm{O_h} \f$.
	QPMS_PGS_I,  ///< Chiral icosahedral symmetry \f$ \mathrm{I} \f$.
	QPMS_PGS_IH,  ///< Full icosahedral symmetry \f$ \mathrm{I_h} \f$.

	// Continuous axial groups
	QPMS_PGS_CINF, ///< \f$ \mathrm{C_\infty} \f$
	QPMS_PGS_CINFH, ///< \f$ \mathrm{C_{\infty h}} \f$
	QPMS_PGS_CINFV, ///< \f$ \mathrm{C_{\infty v}} \f$
	QPMS_PGS_DINF, ///< \f$ \mathrm{D_\infty} \f$
	QPMS_PGS_DINFH, ///< \f$ \mathrm{D_{\infty h}} \f$

	// Remaining continuous groups
	QPMS_PGS_SO3, 	///< Special orthogonal group \f$ \mathrm{SO(3)}.
	QPMS_PGS_O3, 	///< Orthogonal group \f$ \mathrm{O(3)}.
} qpms_pointgroup_class;

/// Full characterisation of a 3D point group.
typedef struct qpms_pointgroup_t {
	qpms_pointgroup_class c; ///< Point group classification.
	qpms_gmi_t n; ///< Order of the rotational subgroup \f$ \mathrm{C_n} \f$ of finite axial groups.
	/// Transformation between this point group and the "canonical" point group of the same type.
	/**
	 * Each 3D point group of a given type (specified by the \a c 
	 * and \a n members) has its isomorphous "canonical instance", 
	 * typically with the main rotation axis identical to the \a z 
	 * cartesian axis and a mirror plane (if applicable) orthogonal 
	 * to the \a x cartesian axis.
	 *
	 * If \f$ o \f$ is a transformation specified by \a orientation, 
	 * then an element \f$ g \f$ of this group can be written
	 * as \f[
	 * 	g = o g_\mathrm{can.} o^{-1}
	 * \f] where \f$ g_\mathrm{can.} \f$ is a corresponding element
	 * from the canonical instance.
	 *
	 * CHECKME \f$ o \f$ transforms the cartesian \a z axis to the
	 * main rotation axis of this group (if applicable).
	 *
	 * TODO more detailed specification about the canonical instances
	 * and in which direction \a orientation goes.
	 */
	qpms_irot3_t orientation;
} qpms_pointgroup_t;


/// An abstract T-matrix without actual elements, but with info about particle symmetry.
typedef struct qpms_abstract_tmatrix_t {
        const qpms_vswf_set_spec_t *spec;
	qpms_pointgroup_t sym;
} qpms_abstract_tmatrix_t;


/// A type holding electric permittivity and magnetic permeability of a material.
typedef struct qpms_epsmu_t {
	complex double eps; ///< Relative permittivity.
	complex double mu; ///< Relative permeability.
} qpms_epsmu_t;

struct qpms_tolerance_spec_t; // See tolerances.h

#define lmcheck(l,m) assert((l) >= 1 && abs(m) <= (l))
#endif // QPMS_TYPES
