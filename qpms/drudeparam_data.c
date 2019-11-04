#include "materials.h"
#include <gsl/gsl_const_mksa.h>

// Data from drudelorentz.com

#define EH (GSL_CONST_MKSA_ELECTRON_VOLT / GSL_CONST_MKSA_PLANCKS_CONSTANT_HBAR)

static const qpms_ldparams_t LDPARAMS_AU = {
  1, // eps_inf
  9.03*EH, // omega_p
  6, // n
  {
    {0.76 ,  0.   *EH, 0.053*EH},
    {0.024,  0.415*EH, 0.241*EH},
    {0.01 ,  0.83 *EH, 0.345*EH},
    {0.071,  2.969*EH, 0.87 *EH},
    {0.601,  4.304*EH, 2.494*EH},
    {4.384, 13.32 *EH, 2.214*EH}
  }
};
const qpms_ldparams_t *const QPMS_LDPARAMS_AU = &LDPARAMS_AU;

static const qpms_ldparams_t LDPARAMS_AG = {
  1, // eps_inf
  9.01*EH, // omega_p
  6, // n
  {
    {0.84 ,  0.   *EH, 0.053*EH},
    {0.065,  0.816*EH, 3.886*EH},
    {0.124,  4.481*EH, 0.452*EH},
    {0.111,  8.185*EH, 0.065*EH},
    {0.84 ,  9.083*EH, 0.916*EH},
    {5.646, 20.29 *EH, 2.419*EH}
  }
};
const qpms_ldparams_t *const QPMS_LDPARAMS_AG = &LDPARAMS_AG;

static const qpms_ldparams_t LDPARAMS_CU = {
  1, // eps_inf
  10.83*EH, // omega_p
  5, // n
  {
    {0.575,  0.   *EH, 0.03 *EH},
    {0.061,  0.291*EH, 0.378*EH},
    {0.104,  2.957*EH, 1.056*EH},
    {0.723,  5.3  *EH, 3.213*EH},
    {0.638, 11.18 *EH, 4.305*EH}
  }
};
const qpms_ldparams_t *const QPMS_LDPARAMS_CU = &LDPARAMS_CU;

static const qpms_ldparams_t LDPARAMS_AL = {
  1, // eps_inf
  14.98*EH, // omega_p
  5, // n
  {
    {0.523, 0.   *EH, 0.047*EH},
    {0.227, 0.162*EH, 0.333*EH},
    {0.05 , 1.544*EH, 0.312*EH},
    {0.166, 1.808*EH, 1.251*EH},
    {0.03 , 3.473*EH, 3.382*EH},
  }
};
const qpms_ldparams_t *const QPMS_LDPARAMS_AL = &LDPARAMS_AL;

static const qpms_ldparams_t LDPARAMS_CR = {
  1, // eps_inf
  10.75*EH, // omega_p
  5, // n
  {
    {0.168, 0.   *EH, 0.047*EH},
    {0.151, 0.121*EH, 3.175*EH},
    {0.15 , 0.543*EH, 1.305*EH},
    {1.149, 1.97 *EH, 2.676*EH},
    {0.825, 8.775*EH, 1.335*EH},
  }
};
const qpms_ldparams_t *const QPMS_LDPARAMS_CR = &LDPARAMS_CR;
	
static const qpms_ldparams_t LDPARAMS_TI = {
  1, // eps_inf
  7.29*EH, // omega_p
  5, // n
  {
    {0.148,  0.   *EH, 0.082*EH},
    {0.899,  0.777*EH, 2.276*EH},
    {0.393,  1.545*EH, 2.518*EH},
    {0.187,  2.509*EH, 1.663*EH},
    {0.001, 19.43 *EH, 1.762*EH},
  }
};
const qpms_ldparams_t *const QPMS_LDPARAMS_TI = &LDPARAMS_TI;
	
static const qpms_ldparams_t LDPARAMS_BE = {
  1, // eps_inf
  18.51*EH, // omega_p
  5, // n
  {
    {0.035, 0.   *EH,	0.035*EH},
    {0.031, 0.1  *EH,	1.664*EH},
    {0.14 , 1.032*EH,	3.395*EH},
    {0.53 , 3.183*EH,	4.454*EH},
    {0.13 , 4.604*EH,	1.802*EH},
  }
};
const qpms_ldparams_t *const QPMS_LDPARAMS_BE = &LDPARAMS_BE;

static const qpms_ldparams_t LDPARAMS_NI = {
  1, // eps_inf
  10.75*EH, // omega_p
  5, // n
  {
    {0.168, 0.   *EH, 0.168*EH},
    {0.151, 0.121*EH, 3.175*EH},
    {0.15 , 0.543*EH, 1.305*EH},
    {1.149, 1.97 *EH, 2.676*EH},
    {0.825, 8.775*EH, 1.335*EH},
  }
};
const qpms_ldparams_t *const QPMS_LDPARAMS_NI = &LDPARAMS_NI;

static const qpms_ldparams_t LDPARAMS_PD = {
  1, // eps_inf
  9.72*EH, // omega_p
  5, // n
  {
    {0.33 , 0.   *EH, 0.008*EH},
    {0.649, 0.336*EH, 2.95 *EH},
    {0.121, 0.501*EH, 0.555*EH},
    {0.638, 1.659*EH, 4.621*EH},
    {0.453, 5.715*EH, 3.236*EH},
  }
};
const qpms_ldparams_t *const QPMS_LDPARAMS_PD = &LDPARAMS_PD;

static const qpms_ldparams_t LDPARAMS_PT = {
  1, // eps_inf
  9.59*EH, // omega_p
  5, // n
  {
    {0.333, 0.   *EH, 0.333*EH},
    {0.191, 0.78 *EH, 0.517*EH},
    {0.659, 1.314*EH, 1.838*EH},
    {0.547, 3.141*EH, 3.668*EH},
    {3.576, 9.249*EH, 8.517*EH},
  }
};
const qpms_ldparams_t *const QPMS_LDPARAMS_PT = &LDPARAMS_PT;

static const qpms_ldparams_t LDPARAMS_W = {
  1, // eps_inf
  13.22*EH, // omega_p
  5, // n
  {
    {0.206, 0.   *EH, 0.206 *EH},
    {0.054, 1.004*EH, 0.53  *EH},
    {0.166, 1.917*EH, 1.281 *EH},
    {0.706, 3.58 *EH, 3.332 *EH},
    {2.59 , 7.498*EH, 5.836 *EH},
  }
};
const qpms_ldparams_t *const QPMS_LDPARAMS_W = &LDPARAMS_W;
