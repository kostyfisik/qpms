#include "lattices.h"
#include <limits.h>
#include <math.h>

// generic converting extractors

PGenPolReturnData PGen_next_pol_from_cart2(PGen *g) {
  const PGenCart2ReturnData c = PGen_next_cart2(g);
  if (c.flags & PGEN_DONE)
    return PGenPolDoneVal;
  else {
    PGenPolReturnData p;
    p.flags = (c.flags & ~QPMS_COORDS_BITRANGE) | QPMS_COORDS_POL;
    p.point_pol = cart2pol(c.point_cart2);
    return p;
  }
}

PGenCart2ReturnData PGen_next_cart2_from_pol(PGen *g) {
  const PGenPolReturnData p = PGen_next_pol(g);
  if (p.flags & PGEN_DONE)
    return PGenCart2DoneVal;
  else {
    PGenCart2ReturnData c;
    c.flags = (p.flags & ~QPMS_COORDS_BITRANGE) | QPMS_COORDS_CART2;
    c.point_cart2 = pol2cart(p.point_pol);
    return c;
  }
}

PGenSphReturnData PGen_next_sph_from_cart3(PGen *g) {
  const PGenCart3ReturnData c = PGen_next_cart3(g);
  if (c.flags & PGEN_DONE)
    return PGenSphDoneVal;
  else {
    PGenSphReturnData s;
    s.flags = (c.flags & ~QPMS_COORDS_BITRANGE) | QPMS_COORDS_SPH;
    s.point_sph = cart2sph(c.point_cart3);
    return s;
  }
}

PGenCart3ReturnData PGen_next_cart3_from_cart2xy(PGen *g) {
  const PGenCart2ReturnData c2 = PGen_next_cart2(g);
  if (c2.flags & PGEN_DONE)
    return PGenCart3DoneVal;
  else {
    PGenCart3ReturnData c3;
    c3.flags = (c2.flags & ~QPMS_COORDS_BITRANGE) | QPMS_COORDS_CART3;
    c3.point_cart3 = cart22cart3xy(c2.point_cart2);
    return c3;
  }
}

PGenSphReturnData PGen_next_sph_from_cart2(PGen *g) {
  const PGenCart2ReturnData c = PGen_next_cart2(g);
  if (c.flags & PGEN_DONE)
    return PGenSphDoneVal;
  else {
    PGenSphReturnData s;
    s.flags = (c.flags & ~QPMS_COORDS_BITRANGE) | QPMS_COORDS_SPH;
    s.point_sph = cart22sph(c.point_cart2);
    return s;
  }
}

PGenCart3ReturnData PGen_next_cart3_from_sph(PGen *g) {
  const PGenSphReturnData s = PGen_next_sph(g);
  if (s.flags & PGEN_DONE)
    return PGenCart3DoneVal;
  else {
    PGenCart3ReturnData c;
    c.flags = (s.flags & ~QPMS_COORDS_BITRANGE) | QPMS_COORDS_CART3;
    c.point_cart3 = sph2cart(s.point_sph);
    return c;
  }
}

// here, various "classes" of the PGenSph point generators are implemented.


// const PGenSphReturnData PGenSphDoneVal = {PGEN_DONE, {0,0,0}}; // defined already in lattices.h
// const PGenCart3ReturnData PGenCart3DoneVal = {PGEN_DONE, {0,0,0}}; // defined already in lattices.h

// General structure of a generator implementation looks like this:

#if 0
//==== PGen_NAME ====

extern const PGenClassInfo PGen_NAME; // forward declaration needed by constructor (may be placed in header file instead)

// Internal state structure
typedef struct PGen_NAME_StateData {
  ...
} PGen_NAME_StateData;

// Constructor
PGenSph PGen_NAME_new(...) {
  g->stateData = malloc(sizeof(PGen_NAME_StateData));
  ...
  PGenSph g = {&PGen_NAME, (void *) stateData};
  return g;
}

// Dectructor
void PGen_NAME_dectructor(PGen *g) {
  ...
  free(g->stateData);
  g->stateData = NULL;
}

// Extractor, spherical coordinate output
PGenSphReturnData PGen_NAME_next_sph(PGen *g) {
  if (g->stateData == NULL) // already destroyed
    return PGenSphDoneVal;
  else {
    PGen_NAME_StateData *s = (PGen_NAME_StateData *) g->stateData;
    if (... /* there are still points to be generated */) {
      ...
      PGenSphReturnData retval = {.../*flags*/, .../*thePoint*/};
      return retval;
    } else {
      PGen_destroy(g);
      return PGenSphDoneVal;
    }
  }
}

// Extractor, 3D cartesian coordinate output
PGenCart3ReturnData PGen_NAME_next_cart3(PGen *g) {
  if (g->stateData == NULL) // already destroyed
    return PGenCart3DoneVal;
  else {
    PGen_NAME_StateData *s = (PGen_NAME_StateData *) g->stateData;
    if (... /* there are still points to be generated */) {
      ...
      PGenCart3ReturnData retval = {.../*flags*/, .../*thePoint*/};
      return retval;
    } else {
      PGen_destroy(g);
      return PGenCart3DoneVal;
    }
  }
}

// Class metadata structure; TODO maybe this can rather be done by macro.
const PGenClassInfo PGen_NAME = {
  "PGen_NAME",
  ?, //dimensionality
  PGEN_COORDS_????, // native coordinate system
  // some of the _next_... fun pointers can be NULL
  PGen_NAME_next,
  PGen_NAME_next_z,
  PGen_NAME_next_pol,
  PGen_NAME_next_sph,
  PGen_NAME_next_cart2,
  PGen_NAME_next_cart3,
  PGen_NAME_fetch,
  PGen_NAME_fetch_z,
  PGen_NAME_fetch_pol,
  PGen_NAME_fetch_sph,
  PGen_NAME_fetch_cart2,
  PGen_NAME_fetch_cart3,
  PGen_NAME_destructor
};

#endif // 0    


//==== PGenSph_FromPoint2DArray ====

// Internal state structure
typedef struct PGen_FromPoint2DArray_StateData {
  const point2d *base;
  size_t len;
  size_t currentIndex;
}PGen_FromPoint2DArray_StateData;

// Constructor
PGen PGen_FromPoint2DArray_new(const point2d *points, size_t len) {
  PGen_FromPoint2DArray_StateData *stateData = malloc(sizeof(PGen_FromPoint2DArray_StateData));
  stateData->base = points;
  stateData->len = len;
  stateData->currentIndex = 0;
  PGen g = {&PGen_FromPoint2DArray, (void *) stateData};
  return g;
}

// Destructor
void PGen_FromPoint2DArray_destructor(PGen *g) {
  free(g->stateData);
  g->stateData = NULL;
}

// Extractor, 2D cartesian (native)
PGenCart2ReturnData PGen_FromPoint2DArray_next_cart2(PGen *g) {
  if (g->stateData == NULL) // already destroyed
    return PGenCart2DoneVal;
  else {
    PGen_FromPoint2DArray_StateData *s = (PGen_FromPoint2DArray_StateData *) g->stateData;
    if (s->currentIndex < s->len) {
      cart2_t thePoint = s->base[s->currentIndex];
      ++(s->currentIndex);
      PGenCart2ReturnData retval = {(PGEN_NOTDONE | PGEN_AT_XY | PGEN_COORDS_CART2), thePoint};
      return retval;
    } else {
      PGen_destroy(g);
      return PGenCart2DoneVal;
    }
  }
}

// Extractor, spherical
PGenSphReturnData PGen_FromPoint2DArray_next_sph(PGen *g) {
  if (g->stateData == NULL) // already destroyed
    return PGenSphDoneVal;
  else {
    PGen_FromPoint2DArray_StateData *s = (PGen_FromPoint2DArray_StateData *) g->stateData;
    if (s->currentIndex < s->len) {
      sph_t thePoint = cart22sph(s->base[s->currentIndex]);
      ++(s->currentIndex);
      PGenSphReturnData retval = {(PGEN_AT_XY | PGEN_COORDS_SPH), thePoint};
      return retval;
    } else {
      PGen_destroy(g);
      return PGenSphDoneVal;
    }
  }
}

const PGenClassInfo PGen_FromPoint2DArray = {
  "PGen_FromPoint2DArray",
  2, // dimensionality
  PGEN_COORDS_CART2,
  NULL,//PGen_FromPoint2DArray_next,
  NULL,
  NULL,//PGen_FromPoint2DArray_next_pol,
  PGen_FromPoint2DArray_next_sph,
  PGen_FromPoint2DArray_next_cart2,
  NULL,//PGen_FromPoint2DArray_next_cart3,
  NULL,//PGen_FromPoint2DArray_fetch,
  NULL,//PGen_FromPoint2DArray_fetch_z,
  NULL,//PGen_FromPoint2DArray_fetch_pol,
  NULL,//PGen_FromPoint2DArray_fetch_sph,
  NULL,//PGen_FromPoint2DArray_fetch_cart2,
  NULL,//PGen_FromPoint2DArray_fetch_cart3,
  PGen_FromPoint2DArray_destructor,
};



//==== PGen_1D ====
//equidistant points along the z-axis;

extern const PGenClassInfo PGen_1D; // forward declaration needed by constructor (may be placed in header file instead)

/* // This had to go to the header file:
enum PGen_1D_incrementDirection{
    //PGEN1D_POSITIVE_INC, // not implemented
    //PGEN1D_NEGATIVE_INC, // not implemented
    PGEN1D_INC_FROM_ORIGIN,
    PGEN1D_INC_TOWARDS_ORIGIN
};
*/

// Internal state structure
typedef struct PGen_1D_StateData {
  long ptindex;
  //long stopindex;
  double minR, maxR;
  bool inc_minR, inc_maxR;
  double a; // lattice period
  double offset; // offset of the zeroth lattice point from origin (will be normalised to interval [-a/2,a/2]
  enum PGen_1D_incrementDirection incdir;
  //bool skip_origin;
} PGen_1D_StateData;

static inline long ptindex_inc(long i) {
  if (i > 0)
    return -i;
  else
    return -i + 1;
}

static inline long ptindex_dec(long i) {
  if (i > 0)
    return -i + 1;
  else
    return -i;
}

// Constructor, specified by maximum and maximum absolute value
PGen PGen_1D_new_minMaxR(double period, double offset, double minR, bool inc_minR, double maxR, bool inc_maxR, 
    PGen_1D_incrementDirection incdir) {
  PGen_1D_StateData *s = malloc(sizeof(PGen_1D_StateData));
  s->minR = minR;
  s->maxR = maxR;
  s->inc_minR = inc_minR;
  s->inc_maxR = inc_maxR;
  s->incdir = incdir;
  period = fabs(period);
  double offset_normalised = offset - period * floor(offset / period); // shift to interval [0, period]
  if (offset_normalised > period / 2) offset_normalised -= period; // and to interval [-period/2, period/2]
  s->offset = offset_normalised;
  if (offset_normalised > 0) // reverse the direction so that the conditions in _next() are hit in correct order
        period *= -1;  
  switch(s->incdir) {
    double curR;
    case PGEN_1D_INC_FROM_ORIGIN:
      s->ptindex = floor(minR / fabs(period));
      while ( (curR = fabs(s->offset + s->ptindex * period)) < minR || (!inc_minR && curR <= minR))
        s->ptindex = ptindex_inc(s->ptindex);
      break;
    case PGEN_1D_INC_TOWARDS_ORIGIN:
      s->ptindex = - ceil(maxR / fabs(period));
      while ( (curR = fabs(s->offset + s->ptindex * period)) > maxR || (!inc_minR && curR >= maxR))
        s->ptindex = ptindex_dec(s->ptindex);
      break;
    default:
      abort(); // invalid argument / not implemented
  }
  s->a = period;
  
  PGen g = {&PGen_1D, (void *) s};
  return g;
}

// Dectructor
void PGen_1D_destructor(PGen *g) {
  free(g->stateData);
  g->stateData = NULL;
}

// Extractor 1D number
PGenZReturnData PGen_1D_next_z(PGen *g) {
  if (g->stateData == NULL) // already destroyed
    return PGenZDoneVal;
  PGen_1D_StateData *s = (PGen_1D_StateData *) g->stateData;
  const double zval = s->ptindex * s->a + s->offset;
  const double r = fabs(zval);
  bool theEnd = false;
  switch (s->incdir) {
    case PGEN_1D_INC_FROM_ORIGIN:
      if (r < s->maxR || (s->inc_maxR && r == s->maxR)) 
        s->ptindex = ptindex_inc(s->ptindex);
      else theEnd = true;
      break; 
    case PGEN_1D_INC_TOWARDS_ORIGIN:
      if (r > s->minR || (s->inc_minR && r == s->minR)) {
        if (s->ptindex == 0) // handle "underflow"
          s->minR = INFINITY;
        else
          s->ptindex = ptindex_dec(s->ptindex);
      } else theEnd = true;
      break;
    default:
        abort(); // invalid value
  }
  if (!theEnd) {
      const PGenZReturnData retval = {PGEN_NOTDONE | PGEN_AT_Z, zval};
      return retval;
  } else {
      PGen_destroy(g);
      return PGenZDoneVal;
  }
}

// Extractor spherical coordinates // TODO remove/simplify
PGenSphReturnData PGen_1D_next_sph(PGen *g) {
  if (g->stateData == NULL) // already destroyed
    return PGenSphDoneVal;
  PGen_1D_StateData *s = (PGen_1D_StateData *) g->stateData;
  const double zval = s->ptindex * s->a + s->offset;
  const double r = fabs(zval);
  bool theEnd = false;
  switch (s->incdir) {
    case PGEN_1D_INC_FROM_ORIGIN:
      if (r < s->maxR || (s->inc_maxR && r == s->maxR)) 
        s->ptindex = ptindex_inc(s->ptindex);
      else theEnd = true;
      break; 
    case PGEN_1D_INC_TOWARDS_ORIGIN:
      if (r > s->minR || (s->inc_minR && r == s->minR)) {
        if (s->ptindex == 0) // handle "underflow"
          s->minR = INFINITY;
        else
          s->ptindex = ptindex_dec(s->ptindex);
      } else theEnd = true;
      break;
    default:
        abort(); // invalid value
  }
  if (!theEnd) {
      const PGenSphReturnData retval = {PGEN_NOTDONE | PGEN_AT_Z | PGEN_COORDS_SPH,
        {r, zval >= 0 ? 0 : M_PI, 0}};
      return retval;
  } else {
      PGen_destroy(g);
      return PGenSphDoneVal;
  }
}

// Class metadata structure; TODO maybe this can rather be done by macro.
const PGenClassInfo PGen_1D = {
  "PGen_1D",
  1, // dimensionality
  PGEN_COORDS_CART1,
  NULL, //PGen_1D_next,
  PGen_1D_next_z,
  NULL,//PGen_1D_next_pol,
  PGen_1D_next_sph,
  NULL,//PGen_1D_next_cart2,
  NULL,//PGen_1D_next_cart3,
  NULL,//PGen_1D_fetch,
  NULL,//PGen_1D_fetch_z,
  NULL,//PGen_1D_fetch_pol,
  NULL,//PGen_1D_fetch_sph,
  NULL,//PGen_1D_fetch_cart2,
  NULL,//PGen_1D_fetch_cart3,
  PGen_1D_destructor
};

//==== PGen_xyWeb ====
// 2D lattice generator in the "spiderweb" style, generated in the "perimetre" order,
// not strictly ordered (or limited) by distance from origin.
// The minR and maxR here refer to the TODO WWHAT

extern const PGenClassInfo PGen_xyWeb; // forward declaration needed by constructor (may be placed in header file instead)

// Internal state structure
typedef struct PGen_xyWeb_StateData {
  long i, j;
  unsigned short phase; // 0 to 5
  long layer;
  long last_layer; // generation stops when layer > last_layer
  double layer_min_height; // this * layer is what minR and maxR are compared to
  double minR, maxR;
  bool inc_minR, inc_maxR;
  cart2_t b1, b2; // lattice vectors
  cart2_t offset; // offset of the zeroth lattice point from origin (TODO will be normalised to the WS cell)
  // TODO type rectangular vs. triangular
  LatticeFlags lf;
} PGen_xyWeb_StateData;

// Constructor
PGen PGen_xyWeb_new(cart2_t b1, cart2_t b2, double rtol, cart2_t offset, double minR, bool inc_minR, double maxR, bool inc_maxR) {
  PGen_xyWeb_StateData *s =  malloc(sizeof(PGen_xyWeb_StateData));
  s->minR = minR; s->maxR = maxR;
  s->inc_minR = inc_minR;
  s->inc_maxR = inc_maxR;
  l2d_reduceBasis(b1, b2, &(s->b1), &(s->b2));
  s->offset = offset; // TODO shorten into the WS cell ?
  s->lf = l2d_detectRightAngles(s->b1, s->b2, rtol);
  s->layer_min_height = l2d_hexWebInCircleRadius(s->b1, s->b2);
  s->layer = ceil(s->minR/s->layer_min_height);
  if(!inc_minR && (s->layer * s->layer_min_height)  <= minR)
    ++(s->layer); 
  s->i = s->layer; s->j = 0; s->phase = 0; // init indices
  s->last_layer = floor(s->maxR/s->layer_min_height);
  if(!inc_maxR && (s->last_layer * s->layer_min_height) >= maxR)
    --(s->last_layer);
  PGen g = {&PGen_xyWeb, (void *) s};
  return g;
}

// Destructor
void PGen_xyWeb_destructor(PGen *g) {
  free(g->stateData);
  g->stateData = NULL;
}

// Extractor (2D cartesian, native)
PGenCart2ReturnData PGen_xyWeb_next_cart2(PGen *g) {
  if (g->stateData == NULL) // already destroyed
    return PGenCart2DoneVal;
  else {
    PGen_xyWeb_StateData * const s = (PGen_xyWeb_StateData *) g->stateData;
    assert(s->layer >= 0);
    if (s->layer <= s->last_layer) {
      const cart2_t thePoint = cart2_add(s->offset,
          cart2_add(cart2_scale(s->i, s->b1), cart2_scale(s->j, s->b2)));
      if(s->layer == 0) { // origin is unique, proceed with next layer
        ++s->layer;
        s->phase = 0;
        s->i = s->layer;
        s->j = 0;
      }
      else if(s->lf & ORTHOGONAL_01) {
        // rectangular or square lattice, four perimeters
        switch(s->phase) {
          case 0: // initial i = l, j = 0
            --s->i;
            ++s->j;
            if(s->i <= 0) ++s->phase;
            break;
          case 1: // initial i = 0, j = l
            --s->i;
            --s->j;
            if(s->j <= 0) ++s->phase;
            break;
          case 2: // initial i = -l, j = 0
            ++s->i;
            --s->j;
            if(s->i >= 0) ++s->phase;
            break;
          case 3: // initial i = 0, j = -l
            ++s->i;
            ++s->j;
            if(s->j >= 0) ++s->phase;
            break;
          default:
            abort();
        }
        if(s->phase == 4) { // phase overflow, start new layer
          ++s->layer;
          s->phase = 0;
          s->i = s->layer;
          s->j = 0;
        }
      } else { // non-rectangular lattice, six perimeters
        switch(s->phase) {
          case 0:
            --s->i;
            ++s->j;
            if(s->i <= 0) ++s->phase;
            break;
          case 1:
            --s->i;
            if(s->i + s->j <= 0) ++s->phase;
            break;
          case 2:
            --s->j;
            if(s->j <= 0) ++s->phase;
            break;
          case 3:
            ++s->i;
            --s->j;
            if(s->i >= 0) ++s->phase;
            break;
          case 4:
            ++s->i;
            if(s->i + s->j >= 0) ++s->phase;
            break;
          case 5:
            ++s->j;
            if(s->j >= 0) ++s->phase;
            break;
          default:
            abort();
        }
        if(s->phase == 6) { // phase overflow, start next layer
          ++s->layer;
          s->phase = 0;
          s->i = s->layer;
          s->j = 0;
        }
      }
      PGenCart2ReturnData retval = {(PGEN_NOTDONE | PGEN_AT_XY | PGEN_COORDS_CART2), thePoint};
      return retval;
    } else {
      PGen_destroy(g);
      return PGenCart2DoneVal;
    }
  }
}

// Class metadata structure; TODO maybe this can rather be done by macro.
const PGenClassInfo PGen_xyWeb = {
  "PGen_xyWeb",
  2,
  PGEN_COORDS_CART2,
  NULL,//PGen_xyWeb_next, // FIXME I should really implement this.
  NULL,//PGen_xyWeb_next_z,
  PGen_next_pol_from_cart2, //NULL,//PGen_xyWeb_next_pol,
  PGen_next_sph_from_cart2, //NULL,//PGen_xyWeb_next_sph,
  PGen_xyWeb_next_cart2, // native
  PGen_next_cart3_from_cart2xy, //NULL,//PGen_xyWeb_next_cart3,
  NULL,//PGen_xyWeb_fetch, // FIXME I should really implement this
  NULL,//PGen_xyWeb_fetch_z,
  NULL,//PGen_xyWeb_fetch_pol,
  NULL,//PGen_xyWeb_fetch_sph,
  NULL,//PGen_xyWeb_fetch_cart2,
  NULL,//PGen_xyWeb_fetch_cart3,
  PGen_xyWeb_destructor
};


size_t PGen_xyWeb_sizecap(cart2_t b1, cart2_t b2, double rtol, cart2_t offset,
    double minR, bool inc_minR, double maxR, bool inc_maxR) 
{
  l2d_reduceBasis(b1, b2, &b1, &b2);
  LatticeFlags lf = l2d_detectRightAngles(b1, b2, rtol);
  double layer_min_height = l2d_hexWebInCircleRadius(b1, b2);
  long layer = ceil(minR / layer_min_height);
  if(!inc_minR && (layer * layer_min_height)  <= minR)
    ++layer; 
  long last_layer = floor(maxR / layer_min_height);
  if(!inc_maxR && (last_layer * layer_min_height) >= maxR)
    --(last_layer);
  // TODO less crude estimate (this one should be safe, however)
  return ((lf & ORTHOGONAL_01) ? 4 : 6) * (last_layer - layer + 1);
}

//==== PGen_LatticeRadialHeap ====

extern const PGenClassInfo PGen_LatticeRadialHeap; // forward declaration needed by constructor (may be placed in header file instead)

// Internal state structure
typedef struct PGen_LatticeRadialHeap_StateData {
  int ldim; // 2 or 3, must 
  int sdim;
  int layer;
  // Minimal distance of a point from the origin in the last layer (if the top of the heap exceeds this, a new layer must be evaluated
  double layer_min_r;
  size_t heap_len;
  size_t heap_capacity;
  double *r_heap;
  int *coord_heap;
  double b[0]; // basis vectors and offset
  double offset_r;
  double minR, maxR;
  bool inc_minR, inc_maxR;
} PGen_LatticeRadialHeap_StateData;

static inline double nd2norm(const double a[], int d) {
  double n = 0;
  for(int i = 0; i < d; ++i) 
    n += a[i]*a[i];
  return sqrt(n);
}

// Constructor
PGen PGen_LatticeRadialHeap_new(int ldim, int sdim, double bvectors[], double offset[], double minR, double maxR,
    bool inc_minR, bool inc_maxR) {
  PGen_LatticeRadialHeap_StateData *s =
       malloc(sizeof(PGen_LatticeRadialHeap_StateData) + (ldim + 1) * sdim * sizeof(double));
  s->ldim = ldim;
  s->sdim = sdim;
  memcpy(s->b, bvectors, ldim * sdim * sizeof(double));
  if (offset) {
    memcpy(s->b + ldim * sdim, offset, sdim * sizeof(double));
    s->offset_r = nd2norm(s->b + ldim*sdim, sdim);
  } else { 
    for (size_t i = 0; i < sdim; ++i) 
      s->b[ldim*sdim + i] = 0;
    s->offset_r = 0;
  }
  s->heap_len = 0;
  s->heap_capacity = 1024;
  QPMS_CRASHING_MALLOC(s->r_heap, sizeof(*s->r_heap) * s->heap_capacity);
  QPMS_CRASHING_MALLOC(s->coord_heap, sizeof(*s->coord_heap) * ldim * s->heap_capacity);
  s->layer = -1;
  s->layer_min_r = -INFINITY;
  s->inc_minR = inc_minR; s->inc_maxR = inc_maxR;
  
  PGen g = {&PGen_LatticeRadialHeap, (void *) s};
  return g;
}

// Increment a counter array with constant sum; in the beginning, both counter[] and counter_cumsum[]
// are expected to be init'd to {thesum, 0, 0, ..., 0}
// Everything is expected to be positive (although the types are signed)
static inline _Bool counter_increment(const int ldim, int counter[], 
    int counter_cumsum[], const int thesum) {
  for(int i = 1; i < ldim; ++i) {
    if(counter_cumsum[i] < thesum) { // Can increment this, do it.
      ++counter[i];
      ++counter_cumsum[i];
      for(int j = i - 1; j > 0; --j) { // Zero the preceding ones
        counter[j] = 0;
        counter_cumsum[j] = counter_cumsum[i];
      }
      // Determine the last digit
      counter[0] = thesum - counter_cumsum[1];
      return 1; // Incremented successfully
    }
  }
  return 0; // Could not increment (counter exhausted)
}

// Assuming the counter is initialised to all non-negative values, step through 
// all the possible sign combinations. In the end, return them back to non-negative
// and return false.
static inline _Bool counter_signcycle(const int ldim, int counter[]) {
  for(int i = 0; i < ldim; ++i) { // Flip signs until a sign flipped to negative (cool, right?)
    counter[i] = -counter[i];
    if (counter[i] < 0)
      return 1;
  }
  return 0;
}

static inline double PGen_LatticeRadialHeap_nextlayer(PGen_LatticeRadialHeap_StateData *s) {
  double mindist = 0;
  s->layer++;
  double minr = +INFINITY;
  const int ldim = s->ldim;
  const int sdim = s->sdim;
  int *counter, *counter_cumsum, *tmp;
  QPMS_CRASHING_CALLOC(counter, s->ldim * 2 * sizeof(*counter));
  counter_cumsum = counter + s->ldim;
  counter[0] = s->layer;
  counter_cumsum[0] = s->layer;

  do {
    if (s->heap_len >= s->heap_capacity - 1) { // Check heap capacity
      s->heap_capacity = s->heap_capacity >= 1024 ? s->heap_capacity * 2 : 1024;
      QPMS_CRASHING_REALLOC(s->r_heap, sizeof(*s->r_heap) * s->heap_capacity);
      QPMS_CRASHING_REALLOC(s->coord_heap, sizeof(*s->coord_heap) * ldim * s->heap_capacity);
    }
    double r = 0;
    for(int i = 0; i < s->sdim; ++i) { // calculate r
      double component = s->b[ldim * sdim + i]; // offset
      for (int j = 0; j < s->ldim; ++j)
        component += s->b[j * sdim + i] * counter[j];
      r += component * component;
    }
    r = sqrt(r);
    minr = MIN(r, minr);
    // Add to the heaps
    int position = s->heap_len++;
    s->r_heap[position] = r;
    memcpy(&s->coord_heap[ldim * position], counter, sizeof(*s->coord_heap) * ldim);
    // bubble-up
    while(position > 0) {
      int parent = (position - 1) / 2;
      if  (s->r_heap[parent] > r)  { // swap
        s->r_heap[position] = s->r_heap[parent];
        memcpy(&s->coord_heap[ldim * position], &s->coord_heap[ldim * parent], sizeof(*s->coord_heap) * ldim);
        s->r_heap[parent] = r;
        memcpy(&s->coord_heap[ldim * parent], counter, sizeof(*s->coord_heap) * ldim);
        position = parent;
      } 
      else break;
    }
  } while (counter_signcycle(s->ldim, counter)
      || counter_increment(s->ldim, counter, counter_cumsum, s->layer));

  free(counter);
  return minr;
}



// sdim-independent generator method
// N.B. This ALWAYS produces, not checking against maxR or destructing the generator itself
// (although it does discard the points with distance smaller (or equal) than minR)
int PGen_LatticeRadialHeap_fillNext(PGen *g, int target[]) {
  if (g->stateData == NULL) // already destroyed
    return -1; //TODO some better error code
  else {
    PGen_LatticeRadialHeap_StateData * const s = (PGen_LatticeRadialHeap_StateData *) g->stateData;
    bool hit = false;
    while(!hit) {
      // Ensure that we have sufficiently filled heap
      while (heap_len < 1 || s->r_heap[0] + s->offset_r > s->layer_min_r) 
        s->layer_min_r = PGen_LatticeRadialHeap_nextlayer(s);
      double r = s->r_heap[0];
      hit = (r > s->minR || (s->inc_minR && r == s->minR));
      if (hit) memcpy(target, coord_heap, s->ldim * sizeof(*target));
      // Heap extract anyway
      // Move last element to root
      --(s->heap_len);
      s->r_heap[0] = s->r_heap[s->heap_len];
      memmove(s->coord_heap, &s->coord_heap[ldim * s->heap_len], ldim * sizeof(*s->coord_heap));
      // Bubble down
      int pos = 0;
      while(1) {
        int largest = pos, kidL = 2*pos+1, kidR = 2*pos+2;
        if (kidL < s->heap_len && s->r_heap[kidL] > s->r_heap[largest])
          largest = kidL;
        if (kidR < s->heap_len && s->r_heap[kidR] > s->r_heap[largest])
          largest = kidR;
        if (largest == pos) 
          break;
        else { // swap
          s->r_heap[pos] = s->r_heap[largest];
          s->r_heap[largest] = r;
          memcpy(&s->coord_heap[ldim * pos], &s->coord_heap[ldim * largest], ldim * sizeof(*s->coord_heap));
          memcpy(&s->coord_heap[ldim * largest], &s->coord_heap[ldim * s->heap_len /*it's still there*/],
              ldim * sizeof(*s->coord_heap));
        }
      }
    }
    return 0;
  }
}


// Destructor
void PGen_LatticeRadialHeap_dectructor(PGen *g) {
  PGen_LatticeRadialHeap_StateData *s = g->stateData;
  free(s->r_heap);
  free(s->coord_heap);
  free(g->stateData);
  g->stateData = NULL;
}

// Extractor, spherical coordinate output
PGenSphReturnData PGen_LatticeRadialHeap_next_sph(PGen *g) {
  if (g->stateData == NULL) // already destroyed
    return PGenSphDoneVal;
  else {
    PGen_LatticeRadialHeap_StateData *s = (PGen_LatticeRadialHeap_StateData *) g->stateData;
    if (... /* there are still points to be generated */) {
      ...
      PGenSphReturnData retval = {.../*flags*/, .../*thePoint*/};
      return retval;
    } else {
      PGen_destroy(g);
      return PGenSphDoneVal;
    }
  }
}

// Extractor, 3D cartesian coordinate output
PGenCart3ReturnData PGen_LatticeRadialHeap_next_cart3(PGen *g) {
  if (g->stateData == NULL) // already destroyed
    return PGenCart3DoneVal;
  else {
    PGen_LatticeRadialHeap_StateData *s = (PGen_LatticeRadialHeap_StateData *) g->stateData;
    if (... /* there are still points to be generated */) {
      ...
      PGenCart3ReturnData retval = {.../*flags*/, .../*thePoint*/};
      return retval;
    } else {
      PGen_destroy(g);
      return PGenCart3DoneVal;
    }
  }
}

// Class metadata structure; TODO maybe this can rather be done by macro.
const PGenClassInfo PGen_LatticeRadialHeap = {
  "PGen_LatticeRadialHeap",
  ?, //dimensionality
  PGEN_COORDS_????, // native coordinate system
  // some of the _next_... fun pointers can be NULL
  PGen_LatticeRadialHeap_next,
  PGen_LatticeRadialHeap_next_z,
  PGen_LatticeRadialHeap_next_pol,
  PGen_LatticeRadialHeap_next_sph,
  PGen_LatticeRadialHeap_next_cart2,
  PGen_LatticeRadialHeap_next_cart3,
  PGen_LatticeRadialHeap_fetch,
  PGen_LatticeRadialHeap_fetch_z,
  PGen_LatticeRadialHeap_fetch_pol,
  PGen_LatticeRadialHeap_fetch_sph,
  PGen_LatticeRadialHeap_fetch_cart2,
  PGen_LatticeRadialHeap_fetch_cart3,
  PGen_LatticeRadialHeap_destructor
};


