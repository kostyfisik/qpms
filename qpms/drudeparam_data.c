#include "materials.h"
#include <gsl/gsl_const_mksa.h>

// Data from drudelorentz.com

#define EH (GSL_CONST_MKSA_ELECTRON_VOLT / GSL_CONST_MKSA_PLANCKS_CONSTANT_HBAR)

static const qpms_ldparams_t LDPARAMS_AU = {
  1, // eps_inf
  9.03*EH, // omega_p
  6, // n
  {
    {0.76 *EH,  0.   *EH, 0.053*EH},
    {0.024*EH,  0.415*EH, 0.241*EH},
    {0.01 *EH,  0.83 *EH, 0.345*EH},
    {0.071*EH,  2.969*EH, 0.87 *EH},
    {0.601*EH,  4.304*EH, 2.494*EH},
    {4.384*EH, 13.32 *EH, 2.214*EH}
  }
};

const qpms_ldparams_t *const QPMS_LDPARAMS_AU = &LDPARAMS_AU;

static const qpms_ldparams_t LDPARAMS_AG = {
  1, // eps_inf
  9.01*EH, // omega_p
  6, // n
  {
    {0.84 *EH,  0.   *EH, 0.053*EH},
    {0.065*EH,  0.816*EH, 3.886*EH},
    {0.124*EH,  4.481*EH, 0.452*EH},
    {0.111*EH,  8.185*EH, 0.065*EH},
    {0.84 *EH,  9.083*EH, 0.916*EH},
    {5.646*EH, 20.29 *EH, 2.419*EH}
  }
};

const qpms_ldparams_t *const QPMS_LDPARAMS_AG = &LDPARAMS_AG;

