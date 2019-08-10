#include <math.h>
#include "qpms_types.h"
#include "qpms_specfunc.h"
#include "gaunt.h"
#include "translations.h"
#include "indexing.h" // TODO replace size_t and int with own index types here
#include <stdbool.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sf_bessel.h>
#include "tiny_inlines.h"
#include "assert_cython_workaround.h"
#include "kahansum.h"
#include <stdlib.h> //abort()
#include <gsl/gsl_sf_coupling.h>
#include "qpms_error.h"
#include "normalisation.h"

//#ifdef QPMS_COMPILE_PYTHON_EXTENSIONS
#include <string.h>

#ifdef QPMS_USE_OMP
#include <omp.h>
#endif

int qpms_cython_trans_calculator_get_AB_arrays_loop(
    const qpms_trans_calculator *c, const qpms_bessel_t J, const int resnd,
    const int daxis, const int saxis,
    char *A_data, const npy_intp *A_shape, const npy_intp *A_strides,
    char *B_data, const npy_intp *B_shape, const npy_intp *B_strides,
    const char *r_data, const npy_intp *r_shape, const npy_intp *r_strides,
    const char *theta_data, const npy_intp *theta_shape, const npy_intp *theta_strides,
    const char *phi_data, const npy_intp *phi_shape, const npy_intp *phi_strides,
    const char *r_ge_d_data, const npy_intp *r_ge_d_shape, const npy_intp *r_ge_d_strides){
  assert(daxis != saxis);
  assert(resnd >= 2);
  int longest_axis = 0;
  int longestshape = 1;
  const npy_intp *resultshape = A_shape, *resultstrides = A_strides;
  // TODO put some restrict's everywhere?
  for (int ax = 0; ax < resnd; ++ax){
    assert(A_shape[ax] == B_shape[ax]);
    assert(A_strides[ax] == B_strides[ax]);
    if (daxis == ax || saxis == ax) continue;
    if (A_shape[ax] > longestshape) {
      longest_axis = ax;
      longestshape = 1;
    }
  }
  const npy_intp longlen = resultshape[longest_axis];

  npy_intp innerloop_shape[resnd];
  for (int ax = 0; ax < resnd; ++ax) {
    innerloop_shape[ax] = resultshape[ax];
  }
  /* longest axis will be iterated in the outer (parallelized) loop. 
   * Therefore, longest axis, together with saxis and daxis, 
   * will not be iterated in the inner loop:
   */  
  innerloop_shape[longest_axis] = 1;
  innerloop_shape[daxis] = 1;
  innerloop_shape[saxis] = 1;

  // these are the 'strides' passed to the qpms_trans_calculator_get_AB_arrays_ext
  // function, which expects 'const double *' strides, not 'char *' ones.
  const npy_intp dstride = resultstrides[daxis] / sizeof(complex double);
  const npy_intp sstride = resultstrides[saxis] / sizeof(complex double);

  int errval = 0;
  // TODO here start parallelisation
  //#pragma omp parallel 
  {
    npy_intp local_indices[resnd];
    memset(local_indices, 0, sizeof(local_indices));
    int errval_local = 0;
    size_t longi;
    //#pragma omp for
    for(longi = 0; longi < longlen; ++longi) {
      // this might be done also in the inverse order, but this is more 
      // 'c-contiguous' way of incrementing the indices
      int ax = resnd - 1;
      while(ax >= 0) {
        /* calculate the correct index/pointer for each array used. 
         * This can be further optimized from O(resnd * total size of 
         * the result array) to O(total size of the result array), but 
         * fick that now
         */
        const char *r_p = r_data + r_strides[longest_axis] * longi;
        const char *theta_p = theta_data + theta_strides[longest_axis] * longi;
        const char *phi_p = phi_data + phi_strides[longest_axis] * longi;
        const char *r_ge_d_p = r_ge_d_data + r_ge_d_strides[longest_axis] * longi;
        char *A_p = A_data + A_strides[longest_axis] * longi;
        char *B_p = B_data + B_strides[longest_axis] * longi;
        for(int i = 0; i < resnd; ++i) {
          // following two lines are probably not needed, as innerloop_shape is there 1 anyway
          // so if i == daxis, saxis, or longest_axis, local_indices[i] is zero.
          if (i == longest_axis) continue;
          if (daxis == i || saxis == i) continue;
          r_p += r_strides[i] * local_indices[i];
          theta_p += theta_strides[i] * local_indices[i];
          phi_p += phi_strides[i] * local_indices[i];
          A_p += A_strides[i] * local_indices[i];
          B_p += B_strides[i] * local_indices[i];
        }

        // perform the actual task here
        errval_local |= qpms_trans_calculator_get_AB_arrays_ext(c, (complex double *)A_p, 
            (complex double *)B_p,
            dstride, sstride,
            // FIXME change all the _ext function types to npy_... so that
            // these casts are not needed
            *((double *) r_p), *((double *) theta_p), *((double *)phi_p),
            (int)(*((npy_bool *) r_ge_d_p)), J);
        if (errval_local) abort();

        // increment the last index 'digit' (ax is now resnd-1; we don't have do-while loop in python)
        ++local_indices[ax];
        while(local_indices[ax] == innerloop_shape[ax] && ax >= 0) {
          // overflow to the next digit but stop when reached below the last one
          local_indices[ax] = 0;
          //local_indices[--ax]++; // dekrementace indexu pod nulu a následná inkrementace poruší paměť FIXME
          ax--;
          if (ax >= 0) local_indices[ax]++;
        }
        if (ax >= 0) // did not overflow, get back to the lowest index
          ax = resnd - 1;
      }
    }
    errval |= errval_local;
  }
  // FIXME when parallelizing
  // TODO Here end parallelisation
  return errval;
}


//#endif // QPMS_COMPILE_PYTHON_EXTENSIONS

