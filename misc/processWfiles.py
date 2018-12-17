#!/usr/bin/env python3

import sys
from qpms import processWfiles_sameKs

npart = int(sys.argv[1])
dest = sys.argv[2]
srcs = sys.argv[3:]

processWfiles_sameKs(srcs, dest, f='d', nparticles=npart)

