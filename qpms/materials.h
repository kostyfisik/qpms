/* \file materials.h
 * \brief Optical properties of materials.
 */
#ifndef QPMS_MATERIALS_H
#define QPMS_MATERIALS_H
#include "qpms_types.h"
#include <gsl/gsl_spline.h>

#ifndef SPEED_OF_LIGHT
/// Speed of light in m/s.
#define SPEED_OF_LIGHT (2.99792458e8)
#endif


/// Gets refractive index of a material from its permeability and permittivity.
/** \f[ n = \sqrt{\mu_r \varepsilon_r} \f] */
static inline complex double qpms_refindex(qpms_epsmu_t em) {
	return csqrt(em.eps * em.mu);
}

/// Gets wave number \a k from angular frequency and material permeability and permittivity.
/** \f[ k = \frac{n\omega}{c_0}  = \frac{\omega\sqrt{\mu_r \varepsilon_r}}{c_0} \f] */
static inline complex double qpms_wavenumber(complex double omega, qpms_epsmu_t em) {
	return qpms_refindex(em)*omega/SPEED_OF_LIGHT;
}

/// Gets (relative) wave impedance \f$ \eta_r \f$ from material permeability and permittivity.
/** \eta_r = \sqrt{\mu_r  / \varepsilon_r} \f] */
static inline complex double qpms_waveimpedance(qpms_epsmu_t em) {
	return csqrt(em.mu / em.eps);
}

/// A \f$ (f_j, \omega_j, \gamma_j) \f for qpms_ldparams_t.
typedef struct qpms_ldparams_triple_t {
	double f;
	double omega;
	double gamma;
} qpms_ldparams_triple_t;

/// Structure holding Lorentz-Drude model parameters of a material.
/** \f[
 * \varepsilon = \varepsilon_\infty + \sum_j=0^{n-1}
 *    \frac{f_j \omega_p^2}{\omega_j^2-\omega^2+i\omega\gamma_j} 
 *  \f]
 */
typedef struct qpms_ldparams_t {
	complex double eps_inf; ///< Permittivity at infinity.
	double omega_p; ///< Plasma frequency.
	size_t n; ///< Number of "oscillators".
	qpms_ldparams_triple_t data[]; ///< "Oscillator" parameters.
} qpms_ldparams_t;

extern const qpms_ldparams_t *const QPMS_LDPARAMS_AG; ///< Lorentz-Drude parameters for silver.
extern const qpms_ldparams_t *const QPMS_LDPARAMS_AU; ///< Lorentz-Drude parameters for gold.

/// Lorentz-Drude permittivity.
complex double qpms_lorentzdrude_eps(const qpms_ldparams_t *, complex double omega);

/// Lorentz-Drude optical properties, with relative permeability set always to one.
qpms_epsmu_t qpms_lorentzdrude_epsmu(const qpms_ldparams_t *, complex double omega);


/// Interpolator of tabulated optical properties.
// TODO use gsl_interp instead of gsl_spline.
typedef struct qpms_permittivity_interpolator_t {
        double *wavelength_m; ///< Wavelength array (in meters).
        double *n; ///< Refraction index array.
        double *k; ///< Attenuation coefficient array.
        gsl_interp *interp_n; ///< Refraction index interpolator object.
        gsl_interp *interp_k; ///< Attenuation coeff interpolator object.
        size_t size; ///< Size of n[], k[], and wavelength_m[].
        // I could add gsl_interp_accel, but that is not necessary.
} qpms_permittivity_interpolator_t;

/// Creates a permittivity interpolator from tabulated wavelengths, refraction indices and extinction coeffs.
qpms_permittivity_interpolator_t *qpms_permittivity_interpolator_create(const size_t incount,
                const double *wavelength_m, ///< Tabulated vacuum wavelength in metres, in strictly increasing order.
                const double *n, ///< Tabulated refraction indices at omega.
                const double *k, ///< Tabulated extinction coefficients.
                const gsl_interp_type *iptype ///< GSL interpolator type
                );

/// Creates a permittivity interpolator from an yml file downloaded from refractiveindex.info website.
qpms_permittivity_interpolator_t *qpms_permittivity_interpolator_from_yml(
                const char *path, ///< Path to the yml file.
                const gsl_interp_type *iptype ///< GSL interpolator type
                );

/// Evaluates interpolated material permittivity at a given angular frequency.
complex double qpms_permittivity_interpolator_eps_at_omega(
                const qpms_permittivity_interpolator_t *interp, double omega_SI);

/// Returns the minimum angular frequency supported by the interpolator.
double qpms_permittivity_interpolator_omega_min(
                const qpms_permittivity_interpolator_t *ip);

/// Returns the minimum angular frequency supported by the interpolator.
double qpms_permittivity_interpolator_omega_max(
                const qpms_permittivity_interpolator_t *ip);

/// Destroy a permittivity interpolator.
void qpms_permittivity_interpolator_free(qpms_permittivity_interpolator_t *interp);

/// Relative permittivity from the Drude model.
static inline complex double qpms_drude_epsilon(
                complex double eps_inf, ///< Relative permittivity "at infinity".
                complex double omega_p, ///< Plasma frequency \f$ \omega_p \f$ of the material.
                complex double gamma_p, ///< Decay constant \f$ \gamma_p \f$ of the material.
                complex double omega ///< Frequency \f$ \omega \f$ at which the permittivity is evaluated.
                ) {
        return eps_inf - omega_p*omega_p/(omega*(omega+I*gamma_p));
}


#endif //QPMS_MATERIALS_H
