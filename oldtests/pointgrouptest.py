#!/usr/bin/env python3
from qpms import PointGroup, CQuat, IRot3, PointGroupClass as PGC

try:
    PointGroup(PGC.CN)
except ValueError:
    print("OK")

try:
    PointGroup(66)
except ValueError:
    print("OK")

C1 = PointGroup(PGC.CN, 1)
print(len(C1))
print(C1.elems())


C2 = PointGroup(PGC.CN, 2)
print(len(C2))
print(C2.elems())

D2H = PointGroup(PGC.DNH, 2)
print(len(D2H))
print(D2H.elems())

exit(0)
