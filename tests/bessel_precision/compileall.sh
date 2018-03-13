#!/bin/bash

c99 -ggdb -o jtest -DJTEST besseltest.c
c99 -ggdb -o ytest -DYTEST besseltest.c
c99 -ggdb -o djtest -DDJTEST besseltest.c
c99 -ggdb -o dytest -DDYTEST besseltest.c
c99 -ggdb -o djtest_steed -DJTEST_STEED besseltest.c

