#include "lattices.h"
#include <limits.h>
#include <math.h>

// here, various "classes" of the PGenSph point generators are implemented.


// const PGenSphReturnData PGenSphDoneVal = {PGEN_DONE, {0,0,0}}; // defined already in lattices.h

// General structure of a generator implementation looks like this:

#if 0
//==== PGenSph_NAME ====

extern const PGenSphClassInfo PGenSph_NAME; // forward declaration needed by constructor (may be placed in header file instead)

// Internal state structure
typedef struct PGenSph_NAME_StateData {
  ...
} PGenSph_NAME_StateData;

// Constructor
PGenSph PGenSph_NAME_new(...) {
  g->stateData = malloc(sizeof(PGenSph_NAME_StateData));
  ...
  PGenSph g = {&PGenSph_NAME, (void *) stateData};
  return g;
}

// Dectructor
void PGenSph_NAME_dectructor(PGenSph *g) {
  ...
  free(g->stateData);
  g->stateData = NULL;
}

// Extractor
PGenSphReturnData PGenSph_NAME_next(PGenSph *g) {
  if (g->stateData == NULL) // already destroyed
    return PGenSphDoneVal;
  else {
    PGenSph_NAME_StateData *s = (PGenSph_NAME_StateData *) g->stateData;
    if (... /* there are still points to be generated */) {
      ...
      PGenSphReturnData retval = {.../*flags*/, .../*thePoint*/};
      return retval;
    } else {
      PGenSph_destroy(g);
      return PGenSphDoneVal;
    }
  }
}

// Class metadata structure; TODO maybe this can rather be done by macro.
const PGenSphClassInfo PGenSph_NAME = {
  "PGenSph_NAME",
  PGenSph_NAME_next,
  PGenSph_NAME_destructor
};

#endif // 0    


//==== PGenSph_FromPoint2DArray ====

// Internal state structure
typedef struct PGenSph_FromPoint2DArray_StateData {
  const point2d *base;
  size_t len;
  size_t currentIndex;
}PGenSph_FromPoint2DArray_StateData;

// Constructor
PGenSph PGenSph_FromPoint2DArray_new(const point2d *points, size_t len) {
  PGenSph_FromPoint2DArray_StateData *stateData = malloc(sizeof(PGenSph_FromPoint2DArray_StateData));
  stateData->base = points;
  stateData->len = len;
  stateData->currentIndex = 0;
  PGenSph g = {&PGenSph_FromPoint2DArray, (void *) stateData};
  return g;
}

// Destructor
void PGenSph_FromPoint2DArray_destructor(PGenSph *g) {
  free(g->stateData);
  g->stateData = NULL;
}

// Extractor
PGenSphReturnData PGenSph_FromPoint2DArray_next(PGenSph *g) {
  if (g->stateData == NULL) // already destroyed
    return PGenSphDoneVal;
  else {
    PGenSph_FromPoint2DArray_StateData *s = (PGenSph_FromPoint2DArray_StateData *) g->stateData;
    if (s->currentIndex < s->len) {
      sph_t thePoint = cart22sph(s->base[s->currentIndex]);
      ++(s->currentIndex);
      PGenSphReturnData retval = {(PGEN_AT_XY | PGEN_NEWR), thePoint};
      return retval;
    } else {
      PGenSph_destroy(g);
      return PGenSphDoneVal;
    }
  }
}

const PGenSphClassInfo PGenSph_FromPoint2DArray = {
  "PGenSph_FromPoint2DArray",
  PGenSph_FromPoint2DArray_next,
  PGenSph_FromPoint2DArray_destructor,
};



//==== PGenSph_zAxis ====
//equidistant points along the z-axis;

extern const PGenSphClassInfo PGenSph_zAxis; // forward declaration needed by constructor (may be placed in header file instead)

/* // This had to go to the header file:
enum PGenSph_zAxis_incrementDirection{
    //PGENSPH_ZAXIS_POSITIVE_INC, // not implemented
    //PGENSPH_ZAXIS_NEGATIVE_INC, // not implemented
    PGENSPH_ZAXIS_INC_FROM_ORIGIN,
    PGENSPH_ZAXIS_INC_TOWARDS_ORIGIN
};
*/

// Internal state structure
typedef struct PGenSph_zAxis_StateData {
  long ptindex;
  //long stopindex;
  double minR, maxR;
  bool inc_minR, inc_maxR;
  double a; // lattice period
  double offset; // offset of the zeroth lattice point from origin (will be normalised to interval [-a/2,a/2]
  enum PGenSph_zAxis_incrementDirection incdir;
  //bool skip_origin;
} PGenSph_zAxis_StateData;

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
PGenSph PGenSph_zAxis_new_minMaxR(double period, double offset, double minR, bool inc_minR, double maxR, bool inc_maxR, 
    PGenSph_zAxis_incrementDirection incdir) {
  PGenSph_zAxis_StateData *s = malloc(sizeof(PGenSph_zAxis_StateData));
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
    case PGENSPH_ZAXIS_INC_FROM_ORIGIN:
      s->ptindex = floor(minR / fabs(period));
      while ( (curR = fabs(s->offset + s->ptindex * period)) < minR || (!inc_minR && curR <= minR))
        s->ptindex = ptindex_inc(s->ptindex);
      break;
    case PGENSPH_ZAXIS_INC_TOWARDS_ORIGIN:
      s->ptindex = - ceil(maxR / fabs(period));
      while ( (curR = fabs(s->offset + s->ptindex * period)) > maxR || (!inc_minR && curR >= maxR))
        s->ptindex = ptindex_dec(s->ptindex);
      break;
    default:
      abort(); // invalid argument / not implemented
  }
  s->a = period;
  
  PGenSph g = {&PGenSph_zAxis, (void *) s};
  return g;
}

// Dectructor
void PGenSph_zAxis_destructor(PGenSph *g) {
  free(g->stateData);
  g->stateData = NULL;
}

// Extractor
PGenSphReturnData PGenSph_zAxis_next(PGenSph *g) {
  if (g->stateData == NULL) // already destroyed
    return PGenSphDoneVal;
  PGenSph_zAxis_StateData *s = (PGenSph_zAxis_StateData *) g->stateData;
  const double zval = s->ptindex * s->a + s->offset;
  const double r = fabs(zval);
  bool theEnd = false;
  switch (s->incdir) {
    case PGENSPH_ZAXIS_INC_FROM_ORIGIN:
      if (r < s->maxR || (s->inc_maxR && r == s->maxR)) 
        s->ptindex = ptindex_inc(s->ptindex);
      else theEnd = true;
      break; 
    case PGENSPH_ZAXIS_INC_TOWARDS_ORIGIN:
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
      const PGenSphReturnData retval = {PGEN_NOTDONE | PGEN_NEWR | PGEN_AT_Z,
        {r, zval >= 0 ? 0 : M_PI, 0}};
      return retval;
  } else {
      PGenSph_destroy(g);
      return PGenSphDoneVal;
  }
}

// Class metadata structure; TODO maybe this can rather be done by macro.
const PGenSphClassInfo PGenSph_zAxis = {
  "PGenSph_zAxis",
  PGenSph_zAxis_next,
  PGenSph_zAxis_destructor
};

#if 0
//==== PGenSph_xyWeb ====
// 2D lattice generator in the "spiderweb" style, generated in the "perimetre" order,
// not strictly ordered (or limited) by distance from origin.
// The minR and maxR here refer to the TODO WWHAT

extern const PGenSphClassInfo PGenSph_xyWeb; // forward declaration needed by constructor (may be placed in header file instead)

// Internal state structure
typedef struct PGenSph_xyWeb_StateData {
  long i, j;
  //long stopindex;
  double minR, maxR;
  bool inc_minR, inc_maxR;
  cart2_t b1, b2; // lattice vectors
  cart2_t offset; // offset of the zeroth lattice point from origin (will be normalised to the WS cell)
} PGenSph_xyWeb_StateData;

// Constructor
PGenSph PGenSph_xyWeb_new(...) {
  g->stateData = malloc(sizeof(PGenSph_xyWeb_StateData));
  ...
  PGenSph g = {&PGenSph_xyWeb, (void *) stateData};
  return g;
}

// Dectructor
void PGenSph_xyWeb_dectructor(PGenSph *g) {
  ...
  free(g->stateData);
  g->stateData = NULL;
}

// Extractor
PGenSphReturnData PGenSph_xyWeb_next(PGenSph *g) {
  if (g->stateData == NULL) // already destroyed
    return PGenSphDoneVal;
  else {
    PGenSph_xyWeb_StateData *s = (PGenSph_xyWeb_StateData *) g->stateData;
    if (... /* there are still points to be generated */) {
      ...
      PGenSphReturnData retval = {.../*flags*/, .../*thePoint*/};
      return retval;
    } else {
      PGenSph_destroy(g);
      return PGenSphDoneVal;
    }
  }
}

// Class metadata structure; TODO maybe this can rather be done by macro.
const PGenSphClassInfo PGenSph_xyWeb = {
  "PGenSph_xyWeb",
  PGenSph_xyWeb_next,
  PGenSph_xyWeb_destructor
};

#endif // 0    



#if 0
//==== PGenSph_xyPlane ==== //TODO

extern const PGenSphClassInfo PGenSph_xyPlane; // forward declaration needed by constructor (may be placed in header file instead)

// Internal state structure
typedef struct PGenSph_xyPlane_StateData {
  long i, j;
  //long stopindex;
  double minR, maxR;
  bool inc_minR, inc_maxR;
  cart2_t b1, b2; // lattice vectors
  cart2_t offset; // offset of the zeroth lattice point from origin (will be normalised to the WS cell)
} PGenSph_xyPlane_StateData;

// Constructor
PGenSph PGenSph_xyPlane_new(...) {
  g->stateData = malloc(sizeof(PGenSph_xyPlane_StateData));
  ...
  PGenSph g = {&PGenSph_xyPlane, (void *) stateData};
  return g;
}

// Dectructor
void PGenSph_xyPlane_dectructor(PGenSph *g) {
  ...
  free(g->stateData);
  g->stateData = NULL;
}

// Extractor
PGenSphReturnData PGenSph_xyPlane_next(PGenSph *g) {
  if (g->stateData == NULL) // already destroyed
    return PGenSphDoneVal;
  else {
    PGenSph_xyPlane_StateData *s = (PGenSph_xyPlane_StateData *) g->stateData;
    if (... /* there are still points to be generated */) {
      ...
      PGenSphReturnData retval = {.../*flags*/, .../*thePoint*/};
      return retval;
    } else {
      PGenSph_destroy(g);
      return PGenSphDoneVal;
    }
  }
}

// Class metadata structure; TODO maybe this can rather be done by macro.
const PGenSphClassInfo PGenSph_xyPlane = {
  "PGenSph_xyPlane",
  PGenSph_xyPlane_next,
  PGenSph_xyPlane_destructor
};

#endif // 0    


