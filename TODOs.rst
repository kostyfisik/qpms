URGENT
======

- Check exact normalisation convention of scuff-tmatrix output.
- Check whether the Condon-Shortley phase affects the form of Wigner matrices.
- The xflip, yflip and possible i-factor problem.


TODO label description (obsolete!)
==================================

LMAXVAR - besides to int, the function should support an iterable for the lMax argument or similar, and
          in such case the output should have the according additional dimension. In the end, I want
          to support systems where different nanoparticles have different lMax.

VECTORIZE - add support to process more than one element at once, in general.

TRUEVECTORIZE - remove inner python loops to speed things up

FEATURE - non-urgent general feature to implement
