import warnings
import argparse
#import sys # for debugging purpose, TODO remove in production
import os # because of path
from .types import TMatrixOp, TMatrixSpec, ParticleSpec, LatticeSpec
import collections 

''' # REMOVE IN PRODUCTION
ParticleSpec = collections.namedtuple('ParticleSpec', ['label', 'position', 'tmatrix_spec']) 
TMatrixOp = collections.namedtuple('TMatrixOp',
                ['optype', 'content'])
TMatrixSpec = collections.namedtuple('TMatrixSpec',
        ['lMax_override', 'tmatrix_path', 'ops'])
'''

__TODOs__ = '''
    - Checking validity of T-matrix ops (the arguments of --tr, --sym or similar) according to what is implemented
      in tmatrices.py.
    - Implement a more user-friendly way to define the lattice base vectors and positions of the particles.
      cf. https://stackoverflow.com/questions/2371436/evaluating-a-mathematical-expression-in-a-string/2371789
    - low priority: allow to perform some more custom operations on T-Matrix, using some kind of parsing from the previous point
    - Autodetect symmetries

'''

def make_action_sharedlist(opname, listname):
    class opAction(argparse.Action):
        def __call__(self, parser, args, values, option_string=None):
            if (not hasattr(args, listname)) or getattr(args, listname) is None:
                setattr(args, listname, list())
            getattr(args,listname).append((opname, values))
    return opAction


def add_argparse_k_output_options(parser):
    parser.add_argument('--kdensity', '--k_density', action='store', type=int, nargs='+', default=33, help='Number of k-points per x-axis segment FIXME DESCRIPTION')
    parser.add_argument('--bz_coverage', action='store', type=float, default=1., help='Brillouin zone coverage in relative length (default 1 for whole 1. BZ)')
    parser.add_argument('--bz_edge_width', action='store', type=float, default=0., help='Width of the more densely covered belt along the 1. BZ edge in relative lengths')
    parser.add_argument('--bz_edge_factor', action='store', type=float, default=8., help='Relative density of the belt along the 1. BZ edge w.r.t. k_density (default==8)')
    parser.add_argument('--bz_edge_twoside', action='store_true', help='Compute also the parts of the densely covered edge belt outside the 1. BZ')
    parser.add_argument('--bz_corner_width', action='store', type=float, default=0., help='Size of the more densely covered subcell along the 1. BZ corners in relative lengths')
    parser.add_argument('--bz_corner_factor', action='store', type=float, default=16., help='Relative density of the subcell along the 1. BZ corner w.r.t. k_density (default==16)')
    parser.add_argument('--bz_corner_twoside', action='store_true', help='Compute also the parts of the densely covered subcell outside the 1. BZ')
    return

def add_argparse_unitcell_definitions(parser):
    #TODO create some user-friendlier way to define lattice vectors, cf. https://stackoverflow.com/questions/2371436/evaluating-a-mathematical-expression-in-a-string/2371789
    parser.add_argument('--lattice_base', nargs=4, action='store', type=float, required=True, help='Lattice basis vectors x1, y1, x2, y2')
    parser.add_argument('--particle', '-p', nargs='+', action=make_action_sharedlist('particle', 'particlespec'), help='Particle label, coordinates x,y, and (optionally) path to the T-Matrix.')
    parser.add_argument('--TMatrix', '-t', nargs='+', action=make_action_sharedlist('TMatrix_path', 'particlespec'), help='Path to TMatrix file')
    parser.add_argument('--background_permittivity', action='store', type=float, default=1., help='Background medium relative permittivity (default 1)')
    parser.add_argument('--lMax', action=make_action_sharedlist('lMax', 'particlespec'), nargs='+', help='Override lMax from the TMatrix file')
    popgrp=parser.add_argument_group(title='Operations')
    popgrp.add_argument('--tr', dest='ops', nargs='+', action=make_action_sharedlist('tr', 'ops'), default=list()) # the default value for dest can be set once
    popgrp.add_argument('--sym', dest='ops', nargs='+', action=make_action_sharedlist('sym', 'ops'))
    #popgrp.add_argument('--mult', dest='ops', nargs=3, metavar=('INCSPEC', 'SCATSPEC', 'MULTIPLIER'), action=make_action_sharedlist('mult', 'ops'))
    #popgrp.add_argument('--multl', dest='ops', nargs=3, metavar=('INCL[,INCL,...]', 'SCATL[,SCATL,...]', 'MULTIPLIER'), action=make_action_sharedlist('multl', 'ops'))
    return

def add_argparse_infinite_lattice_options(parser):
    parser.add_argument('--maxlayer', action='store', type=int, default=100, help='How far to sum the lattice points to obtain the dispersion')
    parser.add_argument('--gaussian', action='store', type=float, metavar='σ', help='Use a gaussian envelope for weighting the interaction matrix contributions (depending on the distance), measured in unit cell lengths (?) FIxME).')
    return

def add_argparse_output_options(parser):
    parser.add_argument('--output_prefix', action='store', required=True, help='Prefix to the npz output (will be appended frequency, hexside and chunkno)')
    parser.add_argument('--plot_TMatrix', action='store_true', help='Visualise TMatrix on the first page of the output')
    parser.add_argument('--scp_to', action='store', metavar='N', type=str, help='SCP the output files to a given destination')
    parser.add_argument('--chunklen', action='store', type=int, default=1000, help='Number of k-points per output file (default 1000)')
    return

def add_argparse_common_options(parser):
    parser.add_argument('--eVfreq', action='store', required=True, type=float, help='Frequency in eV')
    parser.add_argument('--verbose', '-v', action='count', help='Be verbose (about computation times, mostly)')
    parser.add_argument('--frequency_multiplier', action='store', type=float, default=1., help='Multiplies the frequencies in the TMatrix file by a given factor.')

def arg_preprocess_particles(pargs, d=None, return_tuple=False):
    ''' 
    Nanoparticle position and T-matrix path parsing 

    parser: ArgumentParser on which add_argparse_unitcell_definitions() and whose
    parse_args() has been called.

    returns a list of ParticleSpec objects
    '''
    TMatrix_paths = dict()
    lMax_overrides = dict()
    default_TMatrix_path = None
    default_lMax_override = None
    if not any((arg_type == 'particle') for (arg_type, arg_content) in pargs.particlespec):
        # no particles positions given: suppose only one per unit cell, in the cell origin
        positions = {None: (0.0)}
    else:
        positions = dict()
    for arg_type, arg_content in pargs.particlespec:
        if arg_type == 'particle': # --particle option
            if  3 <= len(arg_content) <= 4:
                try:
                    positions[arg_content[0]] = (float(arg_content[1]), float(arg_content[2]))
                except ValueError as e:
                    e.args += ("second and third argument of --particle must be valid floats, given: ", arg_content)
                    raise
                if len(arg_content) == 4:
                    if arg_content[0] in TMatrix_paths:
                        warnings.warn('T-matrix path for particle \'%s\' already specified.' 
                        'Overriding with the last value.' % arg_content[0], SyntaxWarning)
                    TMatrix_paths[arg_content[0]] = arg_content[3]

            else:
                raise ValueError("--particle expects 3 or 4 arguments, %d given: " % len(arg_content), arg_content)
        elif arg_type == 'TMatrix_path': # --TMatrix option
            if len(arg_content) == 1: # --TMatrix default_path
                if default_TMatrix_path is not None:
                    warnings.warn('Default T-matrix path already specified. Overriding with the last value.', SyntaxWarning)
                default_TMatrix_path = arg_content[0]
            elif len(arg_content) > 1: # --TMatrix label [label2 [...]] path
                for label in arg_content[:-1]:
                    if label in TMatrix_paths.keys():
                        warnings.warn('T-matrix path for particle \'%s\' already specified.' 
                        'Overriding with the last value.' % label, SyntaxWarning)
                    TMatrix_paths[label] = arg_content[-1]
        elif arg_type == 'lMax': # --lMax option
            if len(arg_content) == 1: # --lMax default_lmax_override
                if default_lMax_override is not None:
                    warnings.warn('Default lMax override value already specified. Overriding the last value.', SyntaxWarning)
                default_lMax_override = int(arg_content[-1])
            else:
                for label in arg_content[:-1]:
                    if label in lMax_overrides.keys:
                        warnings.warn('lMax override for particle \'%s\' already specified.'
                            'overriding with the last value.' % label, SyntaxWarning)
                    lMax_overrides[label] = int(arg_content[-1])
        else: assert False, 'unknown option type'
    # Check the info from positions and TMatrix_paths and lMax_overrides
    if not set(TMatrix_paths.keys()) <= set(positions.keys()):
        raise ValueError("T-Matrix path(s) for particle(s) labeled %s was given, but not their positions" 
                % str(set(TMatrix_paths.keys()) - set(positions.keys())))
    if not set(lMax_overrides.keys()) <= set(positions.keys()):
        raise ValueError("lMax override(s) for particle(s) labeled %s was given, but not their positions"
                %str(set(lMax_overrides.keys()) - set(positions.keys())))
    if (set(TMatrix_paths.keys()) != set(positions.keys())) and default_TMatrix_path is None:
        raise ValueError("Position(s) of particles(s) labeled %s was given without their T-matrix"
            " and no default T-matrix was specified" 
            % str(set(positions.keys()) - set(TMatrix_paths.keys())))
    # Fill default_TMatrix_path to those that don't have its own
    for label in (set(positions.keys()) - set(TMatrix_paths.keys())):
        TMatrix_paths[label] = default_TMatrix_path
    for path in TMatrix_paths.values():
        if not os.path.exists(path):
            raise ValueError("Cannot access T-matrix file %s. Does it exist?" % path)

    # Assign (pre-parse) the T-matrix operations to individual particles
    ops = dict()
    for label in positions.keys(): ops[label] = list()
    for optype, arg_content in pargs.ops:
        # if, no label given, apply to all, otherwise on the specifield particles
        for label in (positions.keys() if len(arg_content) == 1 else arg_content[:-1]): 
            try:
                ops[label].append(TMatrixOp(optype, arg_content[-1]))
            except KeyError as e:
                e.args += 'Specified operation on undefined particle labeled \'%s\'' % label
                raise

    #### Collect all the info about the particles / their T-matrices into one list ####
    # get rid of the non-unique T-matrix specs (so there is only one instance living
    # of each different TMatrixSpec, possibly with multiple references to it
    TMatrix_specs = dict((spec, spec) 
            for spec in (TMatrixSpec(
                    lMax_overrides[label] if label in lMax_overrides.keys() else None, 
                     TMatrix_paths[label], 
                     tuple(ops[label])) 
                for label in positions.keys())
            )
    # particles_specs contains (label, (xpos, ypos), tmspec per element)
    particles_specs = [ParticleSpec(label, positions[label], 
        TMatrix_specs[(lMax_overrides[label] if label in lMax_overrides.keys() else None, 
                       TMatrix_paths[label], 
                       tuple(ops[label]))]
        ) for label in positions.keys()]

    return particles_specs

'''
import argparse, re, random, string
import subprocess
from scipy.constants import hbar, e as eV, pi, c
import warnings

parser = argparse.ArgumentParser()
pargs=parser.parse_args()
print(pargs)

exit(0) ###

maxlayer=pargs.maxlayer
#DEL hexside=pargs.hexside
eVfreq = pargs.eVfreq
freq = eVfreq*eV/hbar
verbose=pargs.verbose

#DEL TMatrix_file = pargs.TMatrix

epsilon_b = pargs.background_permittivity #2.3104
gaussianSigma = pargs.gaussian if pargs.gaussian else None # hexside * 222 / 7
interpfreqfactor = pargs.frequency_multiplier
scp_dest = pargs.scp_to if pargs.scp_to else None
kdensity = pargs.kdensity
chunklen = pargs.chunklen

# -----------------finished basic CLI parsing (except for op arguments) ------------------
from qpms.timetrack import _time_b, _time_e
btime=_time_b(verbose)

import qpms
import numpy as np
import os, sys, warnings, math
from scipy import interpolate
nx = None
s3 = math.sqrt(3)


# specifikace T-matice zde
cdn = c/ math.sqrt(epsilon_b)
TMatrices_orig, freqs_orig, freqs_weirdunits_orig, lMaxTM = qpms.loadScuffTMatrices(TMatrix_file)
lMax = lMaxTM
if pargs.lMax:
    lMax = pargs.lMax if pargs.lMax else lMaxTM    
my, ny = qpms.get_mn_y(lMax)
nelem = len(my)
if pargs.lMax: #force commandline specified lMax
    TMatrices_orig = TMatrices_orig[...,0:nelem,:,0:nelem]

TMatrices = np.array(np.broadcast_to(TMatrices_orig[:,nx,:,:,:,:],(len(freqs_orig),2,2,nelem,2,nelem)) )

#TMatrices[:,:,:,:,:,ny==3] *= factor13inc
#TMatrices[:,:,:,ny==3,:,:] *= factor13scat
xfl = qpms.xflip_tyty(lMax)
yfl = qpms.yflip_tyty(lMax)
zfl = qpms.zflip_tyty(lMax)
c2rot = qpms.apply_matrix_left(qpms.yflip_yy(3),qpms.xflip_yy(3),-1)

reCN = re.compile('(\d*)C(\d+)')
#TODO C nekonečno

for op in ops:
    if op[0] == 'all':
        targets = (0,1)
    elif isinstance(op[0],int):
        targets = (op[0],)
    else:
        targets = op[0]
        
    if op[1] == 'sym':
        mCN = reCN.match(op[2]) # Fuck van Rossum for not having assignments inside expressions
        if op[2] == 'σ_z':
            for t in targets:
                TMatrices[:,t] = (TMatrices[:,t] + qpms.apply_ndmatrix_left(zfl,qpms.apply_ndmatrix_left(zfl, TMatrices[:,t], (-4,-3)),(-2,-1)))/2
        elif op[2] == 'σ_y':
            for t in targets:
                TMatrices[:,t] = (TMatrices[:,t] + qpms.apply_ndmatrix_left(yfl,qpms.apply_ndmatrix_left(yfl, TMatrices[:,t], (-4,-3)),(-2,-1)))/2
        elif op[2] == 'σ_x':
            for t in targets:
                TMatrices[:,t] = (TMatrices[:,t] + qpms.apply_ndmatrix_left(xfl,qpms.apply_ndmatrix_left(xfl, TMatrices[:,t], (-4,-3)),(-2,-1)))/2
        elif op[2] == 'C2': # special case of the latter
            for t in targets:
                TMatrices[:,t] = (TMatrices[:,t] + qpms.apply_matrix_left(c2rot,qpms.apply_matrix_left(c2rot, TMatrices[:,t], -3),-1))/2                
        elif mCN:
            rotN = int(mCN.group(2))
            TMatrix_contribs = np.empty((rotN,TMatrices.shape[0],2,nelem,2,nelem), dtype=np.complex_)
            for t in targets:
                for i in range(rotN):
                    rotangle = 2*np.pi*i / rotN
                    rot =  qpms.WignerD_yy_fromvector(lMax,np.array([0,0,rotangle]))
                    rotinv = qpms.WignerD_yy_fromvector(lMax,np.array([0,0,-rotangle]))
                    TMatrix_contribs[i] = qpms.apply_matrix_left(rot,qpms.apply_matrix_left(rotinv, TMatrices[:,t], -3),-1)
                TMatrices[:,t] = np.sum(TMatrix_contribs, axis=0) / rotN
        else:
            raise ValueError('\'%d\' is not an implemented symmetry operation' % op[2])
    elif op[1] == 'tr':
        mCN = reCN.match(op[2]) # Fuck van Rossum for not having assignments inside expressions
        if op[2] == 'σ_z':
            for t in targets:
                TMatrices[:,t] = qpms.apply_ndmatrix_left(zfl,qpms.apply_ndmatrix_left(zfl, TMatrices[:,t], (-4,-3)),(-2,-1))
        elif op[2] == 'σ_y':
            for t in targets:
                TMatrices[:,t] = qpms.apply_ndmatrix_left(yfl,qpms.apply_ndmatrix_left(yfl, TMatrices[:,t], (-4,-3)),(-2,-1))
        elif op[2] == 'σ_x':
            for t in targets:
                TMatrices[:,t] = qpms.apply_ndmatrix_left(xfl,qpms.apply_ndmatrix_left(xfl, TMatrices[:,t], (-4,-3)),(-2,-1))
        elif op[2] == 'C2':
            for t in targets:
                TMatrices[:,t] = qpms.apply_matrix_left(c2rot,qpms.apply_matrix_left(c2rot, TMatrices[:,t], -3),-1)       
        elif mCN:
            rotN = int(mCN.group(2))
            power = int(mCN.group(1)) if mCN.group(1) else 1
            TMatrix_contribs = np.empty((rotN,TMatrices.shape[0],2,nelem,2,nelem), dtype=np.complex_)
            for t in targets:
                rotangle = 2*np.pi*power/rotN
                rot = qpms.WignerD_yy_fromvector(lMax, np.array([0,0,rotangle]))
                rotinv = qpms.WignerD_yy_fromvector(lMax, np.array([0,0,-rotangle]))
                TMatrices[:,t] = qpms.apply_matrix_left(rot, qpms.apply_matrix_left(rotinv, TMatrices[:,t], -3),-1)
        else:
            raise ValueError('\'%d\' is not an implemented T-matrix transformation operation' % op[2])
    elif op[1] == 'copy':
        raise # not implemented
    elif op[1] == 'mult':
        raise # not implemented
    elif op[1] == 'multl':
        incy = np.full((nelem,), False, dtype=bool)
        for incl in op[2][0].split(','):
            l = int(incl)
            incy += (l == ny)
        scaty = np.full((nelem,), False, dtype=bool)
        for scatl in op[2][1].split(','):
            l = int(scatl)
            scaty += (l == ny)
        for t in targets:
            TMatrices[np.ix_(np.arange(TMatrices.shape[0]),np.array([t]),np.array([0,1]),scaty,np.array([0,1]),incy)] *= float(op[2][2])
    else:
        raise #unknown operation; should not happen

TMatrices_interp = interpolate.interp1d(freqs_orig*interpfreqfactor, TMatrices, axis=0, kind='linear',fill_value="extrapolate")

klist_full = qpms.generate_trianglepoints(kdensity, v3d=True, include_origin=True)*3*math.pi/(3*kdensity*hexside)
TMatrices_om = TMatrices_interp(freq)

chunkn = math.ceil(klist_full.shape[0] / chunklen)

if verbose:
    print('Evaluating %d k-points in %d chunks' % (klist_full.shape[0], chunkn), file = sys.stderr)
    sys.stderr.flush()

metadata = np.array({
		'lMax' : lMax,
                'maxlayer' : maxlayer,
                'gaussianSigma' : gaussianSigma,
                'epsilon_b' : epsilon_b,
		'hexside' : hexside,
                'chunkn' : chunkn,
                'TMatrix_file' : TMatrix_file,
                'ops' : ops,
                })

for chunki in range(chunkn):
    svdout = '%s_%dnm_%.4f_c%03d.npz' % (pargs.output_prefix, hexside/1e-9, eVfreq, chunki)

    klist = klist_full[chunki * chunklen : (chunki + 1) * chunklen]

    svdres = qpms.hexlattice_zsym_getSVD(lMax=lMax, TMatrices_om=TMatrices_om, epsilon_b=epsilon_b, hexside=hexside, maxlayer=maxlayer,
            omega=freq, klist=klist, gaussianSigma=gaussianSigma, onlyNmin=False, verbose=verbose)

    #((svUfullTElist, svSfullTElist, svVfullTElist), (svUfullTMlist, svSfullTMlist, svVfullTMlist)) = svdres

    np.savez(svdout, omega = freq, klist = klist,
            metadata=metadata,
                     uTE = svdres[0][0],
                     vTE = svdres[0][2],
                     sTE = svdres[0][1],
                     uTM = svdres[1][0],
                     vTM = svdres[1][2],
                     sTM = svdres[1][1],

    )
    svdres=None

    if scp_dest:
        if svdout:
            subprocess.run(['scp', svdout, scp_dest])

_time_e(btime, verbose)
#print(time.strftime("%H.%M:%S",time.gmtime(time.time()-begtime)))
'''
