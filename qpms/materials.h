/* \file materials.h
 * \brief Optical properties of materials.
 */
#ifndef QPMS_MATERIALS_H
#define QPMS_MATERIALS_H
#include "qpms_types.h"
#include <gsl/gsl_spline.h>

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
