#!/usr/bin/env python3
import sys
import argparse

ap = argparse.ArgumentParser()
ap.add_argument("-o", type=str, help='Output file (if not specified, standard output).')
ap.add_argument("input_file", type=str, help="Input file (if not specified, standard input).")

a = ap.parse_args()

if a.o:
    output = open(a.o, 'w')
else:
    output = sys.stdout
if a.input_file:
    input = open(a.input_file, 'r')
else:
    input = sys.stdin

from qpms.constants import eV2SU

for l in input:
    print(eV2SU(float(l)), file=output)

output.close()
input.close()

