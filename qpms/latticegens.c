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
    p.flags = c.flags;
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
    c.flags = p.flags;
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
    s.flags = c.flags;
    s.point_sph = cart2sph(c.point_cart3);
    return s;
  }
}

PGenCart3ReturnData PGen_next_cart3_from_sph(PGen *g) {
  const PGenSphReturnData s = PGen_next_sph(g);
  if (s.flags & PGEN_DONE)
    return PGenCart3DoneVal;
  else {
    PGenCart3ReturnData c;
    c.flags = s.flags;
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
  // some of the _next_... fun pointers can be NULL
  PGen_NAME_next_z,
  PGen_NAME_next_pol,
  PGen_NAME_next_sph,
  PGen_NAME_next_cart2,
  PGen_NAME_next_cart3,
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
      PGenCart2ReturnData retval = {(PGEN_NOTDONE | PGEN_AT_XY | PGEN_NEWR), thePoint};
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
      PGenSphReturnData retval = {(PGEN_AT_XY | PGEN_NEWR), thePoint};
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
  NULL,
  NULL,//PGen_FromPoint2DArray_next_pol,
  PGen_FromPoint2DArray_next_sph,
  PGen_FromPoint2DArray_next_cart2,
  NULL,//PGen_FromPoint2DArray_next_cart3,
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
      const PGenZReturnData retval = {PGEN_NOTDONE | PGEN_NEWR | PGEN_AT_Z,
        zval};
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
      const PGenSphReturnData retval = {PGEN_NOTDONE | PGEN_NEWR | PGEN_AT_Z,
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
  PGen_1D_next_z,
  NULL,//PGen_1D_next_pol,
  PGen_1D_next_sph,
  NULL,//PGen_1D_next_cart2,
  NULL,//PGen_1D_next_cart3,
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
  LatticeType2 lt;
} PGen_xyWeb_StateData;

// Constructor
PGen PGen_xyWeb_new(cart2_t b1, cart2_t b2, double rtol, cart2_t offset, double minR, bool inc_minR, double maxR, bool inc_maxR) {
  PGen_xyWeb_StateData *s =  malloc(sizeof(PGen_xyWeb_StateData));
  s->minR = minR; s->maxR = maxR;
  s->inc_minR = inc_minR;
  s->inc_maxR = inc_maxR;
  l2d_reduceBasis(b1, b2, &(s->b1), &(s->b2));
  s->offset = offset; // TODO shorten into the WS cell ?
  s->lt = l2d_classifyLattice(s->b1, s->b2, rtol);
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
      if(s->layer == 0) // origin is unique, proceed with next layer
        ++s->layer;
      else if(s->lt & (SQUARE | RECTANGULAR)) {
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
      PGenCart2ReturnData retval = {(PGEN_NOTDONE | PGEN_AT_XY | PGEN_NEWR), thePoint};
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
  NULL,//PGen_xyWeb_next_z,
  NULL,//PGen_xyWeb_next_pol,
  NULL,//PGen_xyWeb_next_sph,
  PGen_xyWeb_next_cart2,
  NULL,//PGen_xyWeb_next_cart3,
  PGen_xyWeb_destructor
};


