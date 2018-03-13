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
        print(N(-spherical_bessel_J(l+1,x)
            + l/x*spherical_bessel_J(l,x)), end = ' ', file=file)
    print('', file=file)


