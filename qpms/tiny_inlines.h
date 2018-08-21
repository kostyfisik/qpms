#ifndef TINY_INLINES_H
#define TINY_INLINES_H

static inline int min1pow(int pow) { return (pow % 2) ? -1 : 1; }

// this has shitty precision:
// static inline complex double ipow(int x) { return cpow(I, x); }

static inline complex double ipow(int x) {
  x = ((x % 4) + 4) % 4;
  switch(x) {
    case 0:
      return 1;
    case 1:
      return I;
    case 2:
      return -1;
    case 3:
      return -I;
  }
}

#endif // TINY_INLINES_H
