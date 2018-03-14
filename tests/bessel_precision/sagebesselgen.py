from __future__ import print_function

def printhankelrow(lMax, x, file=sys.stdout):
    print(lMax, N(x), end=' ', file=file);
    for l in range(lMax+1):
        print(N(spherical_hankel1(l,x)), end = ' ', file=file)
    print('', file=file)


def printbesselJrow(lMax, x, file=sys.stdout):
    print(lMax, N(x), end=' ', file=file);
    for l in range(lMax+1):
        print(N(spherical_bessel_J(l,x)), end = ' ', file=file)
    print('', file=file)


def printbesselYrow(lMax, x, file=sys.stdout):
    print(lMax, N(x), end=' ', file=file);
    for l in range(lMax+1):
        print(N(spherical_bessel_Y(l,x)), end = ' ', file=file)
    print('', file=file)

def printbesselDJrow(lMax, x, file=sys.stdout):
    print(lMax, N(x), end=' ', file=file);
    for l in range(lMax+1):
        print(N(-spherical_bessel_J(l+1,x)
            + l/x*spherical_bessel_J(l,x)), end = ' ', file=file)
    print('', file=file)


def printbesselDYrow(lMax, x, file=sys.stdout):
    print(lMax, N(x), end=' ', file=file);
    for l in range(lMax+1):
        print(N(-spherical_bessel_Y(l+1,x)
            + l/x*spherical_bessel_Y(l,x)), end = ' ', file=file)
    print('', file=file)


def genall(lMax):
    f = open('besselDJcases', 'w')
    for o in IntegerRange(1,100):
        printbesselDJrow(lMax, o, file=f)
        printbesselDJrow(lMax, 1/o, file=f)
        printbesselDJrow(lMax, o/sqrt(3), file=f)
    f = open('besselDYcases', 'w')
    for o in IntegerRange(1,100):
        printbesselDYrow(lMax, o, file=f)
        printbesselDYrow(lMax, 1/o, file=f)
        printbesselDYrow(lMax, o/sqrt(3), file=f)
    f = open('besselJcases', 'w')
    for o in IntegerRange(1,100):
        printbesselJrow(lMax, o, file=f)
        printbesselJrow(lMax, 1/o, file=f)
        printbesselJrow(lMax, o/sqrt(3), file=f)
    f = open('besselYcases', 'w')
    for o in IntegerRange(1,100):
        printbesselYrow(lMax, o, file=f)
        printbesselYrow(lMax, 1/o, file=f)
        printbesselYrow(lMax, o/sqrt(3), file=f)

