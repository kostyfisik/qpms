#include <stdio.h>
#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>

#ifndef EPSABS
#define EPSABS 0
#endif
#ifndef EPSREL
#define EPSREL 1e-13
#endif
#ifndef LIMIT
#define LIMIT 30000 //???
#endif
#ifndef R0
#define R0 8e-6
#endif


/* Relevant quadrature methods from gsl:
 * gsl_integration_qagiu
 * ... and that's probably it.
 */

//gsl_function sigma2_integrand;

struct sigma2_integrand_params {
  int n;
  double k, R;
};

static inline double sq(double x) {return x * x;}

double sigma2_integrand(double ksi, void *params) {
  struct sigma2_integrand_params *p = (struct sigma2_integrand_params *) params;
  return exp(-sq(p->R*ksi) + sq(p->k/ksi/2)) * pow(ksi, 2*p->n);
}


int main(int argc, char **argv) {
  struct sigma2_integrand_params p;
  gsl_function F;
  F.function = sigma2_integrand;
  F.params = &p;
  gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(LIMIT);
  gsl_error_handler_t * old_handler=gsl_set_error_handler_off();
  double eta, eta_orig, k_orig, R_orig;
  while (scanf("%d %lf %lf %lf", &(p.n), &k_orig, &R_orig, &eta_orig) == 4) {
    eta = eta_orig * R0;
    p.k = k_orig * R0;
    p.R = R_orig / R0;
    double result, abserr;
    int retval = gsl_integration_qagiu(&F, eta, EPSABS, EPSREL,
        LIMIT, workspace, &result, &abserr);
    double normfac = pow(R0, -2*p.n - 1);
    result *= normfac;
    abserr *= normfac;
    printf("%d %.16g %.16g %.16g %.16g %.16g %d\n", p.n, k_orig, R_orig, 
        eta_orig, result, abserr, retval);
  }
  gsl_integration_workspace_free(workspace);
}


