#include "lattices.h"
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

