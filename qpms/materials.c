#define _POSIX_C_SOURCE 200809L // for getline()
#include <stdio.h>
#include <stddef.h>
#include <unistd.h>
#include <string.h>
#include "materials.h"
#include "qpms_error.h"

#define SQ(x) ((x)*(x))

qpms_permittivity_interpolator_t *qpms_permittivity_interpolator_create(
    const size_t incount, const double *wavelen_m, const double *n, const double *k,
    const gsl_interp_type *iptype)
{
  if (incount <= 0) return NULL;
  qpms_permittivity_interpolator_t *ip;
  QPMS_CRASHING_MALLOC(ip, sizeof(qpms_permittivity_interpolator_t));
  ip->size = incount;
  QPMS_CRASHING_MALLOC(ip->wavelength_m, incount * sizeof(double));
  QPMS_CRASHING_MALLOC(ip->n, incount * sizeof(double));
  QPMS_CRASHING_MALLOC(ip->k, incount * sizeof(double));
  memcpy(ip->wavelength_m, wavelen_m, incount*sizeof(double));
  memcpy(ip->k, k, incount*sizeof(double));
  memcpy(ip->n, n, incount*sizeof(double));
  ip->interp_n = gsl_interp_alloc(iptype, incount);
  ip->interp_k = gsl_interp_alloc(iptype, incount);
  QPMS_ENSURE_SUCCESS(gsl_interp_init(ip->interp_n, ip->wavelength_m, ip->n, incount));
  QPMS_ENSURE_SUCCESS(gsl_interp_init(ip->interp_k, ip->wavelength_m, ip->k, incount));
  return ip;
}

void qpms_permittivity_interpolator_free(qpms_permittivity_interpolator_t *interp)
{
  if(interp) {
    gsl_interp_free(interp->interp_n);
    gsl_interp_free(interp->interp_k);
    free(interp->n);
    free(interp->k);
    free(interp->wavelength_m);
  }
  free(interp);
}

qpms_errno_t qpms_read_refractiveindex_yml(
    FILE *f, ///< file handle
    size_t *const count, ///< Number of successfully loaded triples.
    double* *const lambdas_m, ///< Vacuum wavelengths in metres.
    double* *const n, ///< Read refraction indices.
    double* *const k ///< Read attenuation coeffs.
    )
{
  QPMS_ENSURE(f && lambdas_m && n && k,"f, lambdas_m, n, k are mandatory arguments and must not be NULL.");
  int count_alloc = 128; // First chunk to allocate
  *count = 0;
  QPMS_CRASHING_MALLOC(*lambdas_m, count_alloc * sizeof(double));
  QPMS_CRASHING_MALLOC(*n, count_alloc * sizeof(double));
  QPMS_CRASHING_MALLOC(*k, count_alloc * sizeof(double));
  size_t linebufsz = 256;
  char *linebuf;
  QPMS_CRASHING_MALLOC(linebuf, linebufsz);
  ssize_t readchars;
  bool data_started = false;
  while((readchars = getline(&linebuf, &linebufsz, f)) != -1) {
    if (linebuf[0] == '#') continue;
    // We need to find the beginning of the tabulated data; everything before that is ignored.
    if (!data_started) {
      char *test = strstr(linebuf, "data: |");
      if(test) data_started = true;
      continue;
    }

    if (3 == sscanf(linebuf, "%lf %lf %lf", *lambdas_m + *count, *n + *count , *k + *count)) {
      (*lambdas_m)[*count] *= 1e-6; // The original data is in micrometres.
      ++*count;
      if (*count > count_alloc) {
        count_alloc *= 2;
        QPMS_CRASHING_REALLOC(*lambdas_m, count_alloc * sizeof(double));
        QPMS_CRASHING_REALLOC(*n, count_alloc * sizeof(double));
        QPMS_CRASHING_REALLOC(*k, count_alloc * sizeof(double));
      }
    } else break;
  }
  QPMS_ENSURE(*count > 0, "Could not read any refractive index data; the format must be wrong!");
  free(linebuf);
  return QPMS_SUCCESS;
}

qpms_errno_t qpms_load_refractiveindex_yml(
    const char *path,
    size_t *const count, ///< Number of successfully loaded triples.
    double* *const lambdas_m, ///< Vacuum wavelengths in metres.
    double* *const n, ///< Read refraction indices.
    double* *const k ///< Read attenuation coeffs.
    )
{
  FILE *f = fopen(path, "r");
  QPMS_ENSURE(f, "Could not open refractive index file %s", path);
  qpms_errno_t retval =
    qpms_read_refractiveindex_yml(f, count, lambdas_m, n, k);
  QPMS_ENSURE_SUCCESS(fclose(f));
  return retval;
}

qpms_permittivity_interpolator_t *qpms_permittivity_interpolator_from_yml(
    const char *path, ///< Path to the yml file.
    const gsl_interp_type *iptype ///< GSL interpolator type
    )
{
  size_t count;
  double *lambdas_m, *n, *k;
  QPMS_ENSURE_SUCCESS(qpms_load_refractiveindex_yml(path, &count, &lambdas_m, &n, &k));
  qpms_permittivity_interpolator_t *ip = qpms_permittivity_interpolator_create(
      count, lambdas_m, n, k, iptype);
  free(lambdas_m);
  free(n);
  free(k);
  return ip;
}

complex double qpms_permittivity_interpolator_eps_at_omega(
    const qpms_permittivity_interpolator_t *ip,  double omega_SI)
{
  double lambda, n, k;
  lambda = 2*M_PI*SPEED_OF_LIGHT/omega_SI;
  n = gsl_interp_eval(ip->interp_n, ip->wavelength_m, ip->n, lambda, NULL);
  k = gsl_interp_eval(ip->interp_k, ip->wavelength_m, ip->k, lambda, NULL);
  complex double epsilon = n*n - k*k + 2*n*k*I;
  return epsilon;
}

qpms_epsmu_t qpms_permittivity_interpolator_epsmu_g(
    complex double omega, const void *p) 
{
  const qpms_permittivity_interpolator_t *interp = p;
  static bool imag_already_bitched = false;
  if(cimag(omega) && !imag_already_bitched) 
    QPMS_WARN("Complex frequencies not supported by qpms_permittivity_interpolator_t. Imaginary parts will be discarded!");
  qpms_epsmu_t em;
  em.eps = qpms_permittivity_interpolator_eps_at_omega(interp, omega);
  em.mu = 1;
  return em;
}

double qpms_permittivity_interpolator_omega_max(
    const qpms_permittivity_interpolator_t *ip)
{
  return 2*M_PI*SPEED_OF_LIGHT / ip->wavelength_m[0];
}

double qpms_permittivity_interpolator_omega_min(
    const qpms_permittivity_interpolator_t *ip)
{
  return 2*M_PI*SPEED_OF_LIGHT / ip->wavelength_m[ip->size-1];
}

complex double qpms_lorentzdrude_eps(complex double omega, const qpms_ldparams_t *p) 
{
  complex double eps = 0;
  for(size_t j = 0; j < p->n; ++j) {
    const qpms_ldparams_triple_t d = p->data[j];
    eps += d.f * SQ(p->omega_p) / (SQ(d.omega) - SQ(omega) + I*omega*d.gamma );
  }
  return eps;
}

qpms_epsmu_t qpms_lorentzdrude_epsmu(complex double omega, const qpms_ldparams_t *p) 
{
  qpms_epsmu_t em;
  em.eps = qpms_lorentzdrude_eps(omega, p);
  em.mu = 1;
  return em;
}

qpms_epsmu_t qpms_lorentzdrude_epsmu_g(complex double omega, const void *p) 
{
  return qpms_lorentzdrude_epsmu(omega, (const qpms_ldparams_t *)p);
}

qpms_epsmu_t qpms_epsmu_const_g(complex double omega, const void *p)
{
  return *(const qpms_epsmu_t *)p;
}
