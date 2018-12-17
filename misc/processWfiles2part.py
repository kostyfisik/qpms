#!/usr/bin/env python3

import sys
from qpms import processWfiles_sameKs

dest = sys.argv[1]
srcs = sys.argv[2:]

processWfiles_sameKs(srcs, dest, f='d')

