#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "kahansum.h"
#include "vectors.h"

static const double s3 = 1.732050807568877293527446341505872366942805253810380628055;

// For sorting the points by distance from origin / radius
int cart2_cmpr (const void *p1, const void *p2) {
  const  cart2_t *p1t = (const  cart2_t *)p1;
  const  cart2_t *p2t = (const  cart2_t *)p2;
  double r21 = cart2norm(*p1t);
  double r22 = cart2norm(*p2t);
  if (r21 < r22) return -1;
  else if (r21 > r22) return 1;
  else return 0;
}

typedef struct {
  int npoints; // number of lattice points.
  int capacity; // for how much points memory is allocated
  double maxR; // circle radius, points <= R
   cart2_t *points;
} latticepoints_circle_t;

void sort_cart2points_by_r(cart2_t *points, size_t nmemb) {
  qsort(points, nmemb, sizeof(cart2_t), cart2_cmpr);
}

void latticepoints_circle_free(latticepoints_circle_t *c) {
  free(c->points);
  c->capacity = 0;
}

// "horizontal" orientation of the adjacent A, B points
latticepoints_circle_t generate_hexpoints_hor(double h, double R, cart2_t offset /* if zero, an A is in the origin */) {
  latticepoints_circle_t lat;
  lat.maxR = R;
  lat.npoints = 0;
  int nmax = R / (1.5 * h) + 2; // max no of lattice shifts in each direction (with reserve)
  double unitcellS = s3 * 3 / 2 * h * h; // unit cell area
  lat.capacity = 5 + 2 * (R + 1) * (R + 1) * M_PI / unitcellS; // should be enough with some random reserve
  lat.points = malloc(lat.capacity *sizeof(cart2_t));

  cart2_t BAoffset = {h, 0};

  cart2_t a1 = {-1.5*h, s3/2 *h};
  cart2_t a2 = {1.5*h, s3/2 *h};

  for (int i1 = -nmax; i1 <= nmax; ++i1)
    for (int i2 = -nmax; i2 <= nmax; ++i2) {
      cart2_t Apoint = cart2_add(offset, cart2_add(cart2_scale(i1, a1), cart2_scale(i2, a2)));
      if (cart2norm(Apoint) <= R) {
        assert(lat.npoints < lat.capacity);
        lat.points[lat.npoints] = Apoint;
        lat.npoints++;
      }
      cart2_t Bpoint = cart2_add(Apoint, BAoffset);
      if (cart2norm(Bpoint) <= R) {
        assert(lat.npoints < lat.capacity);
        lat.points[lat.npoints] = Bpoint;
        lat.npoints++;
      }
    }
  sort_cart2points_by_r(lat.points, lat.npoints);
  return lat;
}

latticepoints_circle_t generate_tripoints_ver(double a, double R, cart2_t offset /* if zero, an A is in the origin */) {
  double h = a / s3;
  latticepoints_circle_t lat;
  lat.maxR = R;
  lat.npoints = 0;
  int nmax = R / (1.5 * h) + 2; // max no of lattice shifts in each direction (with reserve)
  double unitcellS = s3 * 3 / 2 * h * h; // unit cell area
  lat.capacity = 5 + (R + 1) * (R + 1) * M_PI / unitcellS; // should be enough with some random reserve
  lat.points = malloc(lat.capacity *sizeof(cart2_t));

  cart2_t a1 = {-1.5*h, s3/2 *h};
  cart2_t a2 = {1.5*h, s3/2 *h};

  for (int i1 = -nmax; i1 <= nmax; ++i1)
    for (int i2 = -nmax; i2 <= nmax; ++i2) {
      cart2_t Apoint = cart2_add(offset, cart2_add(cart2_scale(i1, a1), cart2_scale(i2, a2)));
      if (cart2norm(Apoint) <= R) {
        assert(lat.npoints < lat.capacity);
        lat.points[lat.npoints] = Apoint;
        lat.npoints++;
      }
    }
  sort_cart2points_by_r(lat.points, lat.npoints);
  return lat;
}

latticepoints_circle_t generate_tripoints_hor(double a, double R, cart2_t offset /* if zero, an A is in the origin */) {
  double h = a / s3;
  latticepoints_circle_t lat;
  lat.maxR = R;
  lat.npoints = 0;
  int nmax = R / (1.5 * h) + 2; // max no of lattice shifts in each direction (with reserve)
  double unitcellS = s3 * 3 / 2 * h * h; // unit cell area
  lat.capacity = 5 + (R + 1) * (R + 1) * M_PI / unitcellS; // should be enough with some random reserve
  lat.points = malloc(lat.capacity *sizeof(cart2_t));

  cart2_t a1 = {s3/2 *h, -1.5*h};
  cart2_t a2 = {s3/2 *h, 1.5 * h};

  for (int i1 = -nmax; i1 <= nmax; ++i1)
    for (int i2 = -nmax; i2 <= nmax; ++i2) {
      cart2_t Apoint = cart2_add(offset, cart2_add(cart2_scale(i1, a1), cart2_scale(i2, a2)));
      if (cart2norm(Apoint) <= R) {
        assert(lat.npoints < lat.capacity);
        lat.points[lat.npoints] = Apoint;
        lat.npoints++;
      }
    }
  sort_cart2points_by_r(lat.points, lat.npoints);
  return lat;
}

int main (int argc, char **argv) {
  double h = 1.2;
  cart2_t offset = {0,0};

  latticepoints_circle_t lat = generate_tripoints_ver(h, 30, offset);
  for (int i = 0; i < lat.npoints; ++i) 
    printf("%g %g %g\n", lat.points[i].x, lat.points[i].y, cart2norm(lat.points[i]));
  latticepoints_circle_free(&lat);
}
