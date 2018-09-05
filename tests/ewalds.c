// implementation of the [LT(4.16)] test

#include <qpms/ewald.h>

typedef struct ewaldtest_hex_params {
  qpms_l_t lMax;
  poinnt2d beta; 
  double k; 
  double h;
  double eta;
} ewaldtest_hex_params;


ewaldtest_hex_paraps paramslist = {
  { 3, {1.1, 0.23}, 2.3, 0.97, 0.3},
// end:
  { 0,  {0, 0}, 0, 0, 0}
}



