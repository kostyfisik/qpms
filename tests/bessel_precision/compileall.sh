#!/bin/bash

c99 -ggdb -o jtest -DJTEST besseltest.c -lgsl -lblas -lm
c99 -ggdb -o ytest -DYTEST besseltest.c -lgsl -lblas -lm
c99 -ggdb -o jtest_steed -DJTEST_STEED besseltest.c -lgsl -lblas -lm
c99 -ggdb -o djtest -DDJTEST besseltest.c -lgsl -lblas -lm
c99 -ggdb -o dytest -DDYTEST besseltest.c -lgsl -lblas -lm
c99 -ggdb -o djtest_steed -DDJTEST_STEED besseltest.c -lgsl -lblas -lm

