#include "qpms_specfunc.h"
#include "qpms_types.h"
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_math.h>
#include <stdlib.h>
#include "indexing.h"
#include <string.h>

// Legendre functions also for negative m, see DLMF 14.9.3
qpms_errno_t qpms_legendre_deriv_y_fill(double *target, double *target_deriv, double x, qpms_l_t lMax,
                gsl_sf_legendre_t lnorm, double csphase)
{
        size_t n = gsl_sf_legendre_array_n(lMax);
        	double *legendre_tmp = malloc(n * sizeof(double));
        double *legendre_deriv_tmp = malloc(n * sizeof(double));
        int gsl_errno = gsl_sf_legendre_deriv_array_e(
                        lnorm, (size_t)lMax, x, csphase, legendre_tmp,legendre_deriv_tmp);
        for (qpms_l_t l = 1; l <= lMax; ++l)
                        for (qpms_m_t m = 0; m <= l; ++m) {
                                qpms_y_t y = qpms_mn2y(m,l);
                                size_t i = gsl_sf_legendre_array_index(l,m);
                                target[y] = legendre_tmp[i];
                                target_deriv[y] = legendre_deriv_tmp[i];
                        }
        switch(lnorm) {
                case GSL_SF_LEGENDRE_NONE:
                        for (qpms_l_t l = 1; l <= lMax; ++l)
                                for (qpms_m_t m = 1; m <= l; ++m) {
                                qpms_y_t y = qpms_mn2y(-m,l);
                                size_t i = gsl_sf_legendre_array_index(l,m);
                                // viz DLMF 14.9.3, čert ví, jak je to s cs fasí.
                                double factor = exp(lgamma(l-m+1)-lgamma(l+m+1))*((m%2)?-1:1);
                                target[y] = factor * legendre_tmp[i];
                                target_deriv[y] = factor * legendre_deriv_tmp[i];
                                }
                        break;
                case GSL_SF_LEGENDRE_SCHMIDT:
                case GSL_SF_LEGENDRE_SPHARM:
                case GSL_SF_LEGENDRE_FULL:
                        for (qpms_l_t l = 1; l <= lMax; ++l)
                                for (qpms_m_t m = 1; m <= l; ++m) {
                                qpms_y_t y = qpms_mn2y(-m,l);
                                size_t i = gsl_sf_legendre_array_index(l,m);
                                // viz DLMF 14.9.3, čert ví, jak je to s cs fasí.
                                double factor = ((m%2)?-1:1); // this is the difference from the unnormalised case
                                target[y] = factor * legendre_tmp[i];
                                target_deriv[y] = factor * legendre_deriv_tmp[i];
                                }
                        break;
                default:
                        abort(); //NI
                        break;
        }
        free(legendre_tmp);
        free(legendre_deriv_tmp);
        return QPMS_SUCCESS;
}

qpms_errno_t qpms_legendre_deriv_y_get(double **target, double **dtarget, double x, qpms_l_t lMax, gsl_sf_legendre_t lnorm,
                double csphase)
{

        *target = malloc(sizeof(double)*qpms_lMax2nelem(lMax));
        *dtarget = malloc(sizeof(double)*qpms_lMax2nelem(lMax));
        return qpms_legendre_deriv_y_fill(*target, *dtarget, x, lMax, lnorm, csphase);
}

qpms_pitau_t qpms_pitau_get(double theta, qpms_l_t lMax, qpms_normalisation_t norm)
{
        const double csphase = qpms_normalisation_t_csphase(norm);
        norm = qpms_normalisation_t_normonly(norm);
        qpms_pitau_t res;
        qpms_y_t nelem = qpms_lMax2nelem(lMax);
        res.pi = malloc(nelem * sizeof(double));
        res.tau = malloc(nelem * sizeof(double));
        double ct = cos(theta), st = sin(theta);
        if (1 == fabs(ct)) { // singular case, use DLMF 14.8.2
                memset(res.pi, 0, nelem*sizeof(double));
                memset(res.tau, 0, nelem*sizeof(double));
                res.leg = calloc(nelem, sizeof(double));
                switch(norm) {
                        case QPMS_NORMALISATION_XU:
                                for (qpms_l_t l = 1; l <= lMax; ++l) {
                                        res.leg[qpms_mn2y(0, l)] = (l%2)?ct:1.;
                                        double p = l*(l+1)/2;
                                        const double n = 0.5;
                                        int lpar = (l%2)?-1:1;
                                        res.pi [qpms_mn2y(+1, l)] = -((ct>0) ? -1 : lpar) * p * csphase;
                                        res.pi [qpms_mn2y(-1, l)] = -((ct>0) ? -1 : lpar) * n * csphase;
                                        res.tau[qpms_mn2y(+1, l)] = ((ct>0) ? +1 : lpar) * p * csphase;
                                        res.tau[qpms_mn2y(-1, l)] = -((ct>0) ? +1 : lpar) * n * csphase;
                                }
                                break;
                        case QPMS_NORMALISATION_TAYLOR:
                                for (qpms_l_t l = 1; l <= lMax; ++l) {
                                        res.leg[qpms_mn2y(0, l)] = ((l%2)?ct:1.)*sqrt((2*l+1)*0.25*M_1_PI);
                                        int lpar = (l%2)?-1:1;
                                        double fl = 0.25 * sqrt((2*l+1)*l*(l+1)*M_1_PI);
                                        res.pi [qpms_mn2y(+1, l)] = -((ct>0) ? -1 : lpar) * fl * csphase;
                                        res.pi [qpms_mn2y(-1, l)] = -((ct>0) ? -1 : lpar) * fl * csphase;
                                        res.tau[qpms_mn2y(+1, l)] = ((ct>0) ? +1 : lpar) * fl * csphase;
                                        res.tau[qpms_mn2y(-1, l)] = -((ct>0) ? +1 : lpar) * fl * csphase;
                                }
                                break;
                        case QPMS_NORMALISATION_POWER:
                                for (qpms_l_t l = 1; l <= lMax; ++l) {
                                        res.leg[qpms_mn2y(0, l)] = ((l%2)?ct:1.)*sqrt((2*l+1)/(4*M_PI *l*(l+1)));
                                        int lpar = (l%2)?-1:1;
                                        double fl = 0.25 * sqrt((2*l+1)*M_1_PI);
                                        res.pi [qpms_mn2y(+1, l)] = -((ct>0) ? -1 : lpar) * fl * csphase;
                                        res.pi [qpms_mn2y(-1, l)] = -((ct>0) ? -1 : lpar) * fl * csphase;
                                        res.tau[qpms_mn2y(+1, l)] = ((ct>0) ? +1 : lpar) * fl * csphase;
                                        res.tau[qpms_mn2y(-1, l)] = -((ct>0) ? +1 : lpar) * fl * csphase;

                                }
                                break;
                        default:
                                abort();
                }
        }
        else { // cos(theta) in (-1,1), use normal calculation
                double *legder = malloc(sizeof(double)*qpms_lMax2nelem(lMax));
                res.leg = malloc(sizeof(double)*qpms_lMax2nelem(lMax));
                if (qpms_legendre_deriv_y_fill(res.leg, legder, ct, lMax,
                                        norm == QPMS_NORMALISATION_XU ? GSL_SF_LEGENDRE_NONE
                                                : GSL_SF_LEGENDRE_SPHARM, csphase))
                        abort();
                if (norm == QPMS_NORMALISATION_POWER)
                        /* for Xu (=non-normalized) and Taylor (=sph. harm. normalized)
                         * the correct normalisation is already obtained from gsl_sf_legendre_deriv_array_e().
                         * However, Kristensson ("power") normalisation differs from Taylor
                         * by 1/sqrt(l*(l+1)) factor.
                         */
                        for (qpms_l_t l = 1; l <= lMax; ++l) {
                                double prefac = 1./sqrt(l*(l+1));
                                for (qpms_m_t m = -l; m <= l; ++m) {
                                        res.leg[qpms_mn2y(m,l)] *= prefac;
                                        legder[qpms_mn2y(m,l)] *= prefac;
                                }
                        }
                for (qpms_l_t l = 1; l <= lMax; ++l) {
                        for (qpms_m_t m = -l; m <= l; ++m) {
                                                res.pi [qpms_mn2y(m,l)] = m / st * res.leg[qpms_mn2y(m,l)];
                                                res.tau[qpms_mn2y(m,l)] = - st * legder[qpms_mn2y(m,l)];
                        }
                }
                free(legder);
        }
        res.lMax = lMax;
        return res;
}

void qpms_pitau_free(qpms_pitau_t x) {
        free(x.leg);
        free(x.pi);
        free(x.tau);
}

