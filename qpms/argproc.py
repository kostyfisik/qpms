'''
Common snippets for argument processing in command line scripts; legacy scripts use scripts_common.py instead.
'''

import argparse
import sys
import warnings

def flatten(S):
    if S == []:
        return S
    if isinstance(S[0], list):
        return flatten(S[0]) + flatten(S[1:])
    return S[:1] + flatten(S[1:])

def make_action_sharedlist(opname, listname):
    class opAction(argparse.Action):
        def __call__(self, parser, args, values, option_string=None):
            if (not hasattr(args, listname)) or getattr(args, listname) is None:
                setattr(args, listname, list())
            getattr(args, listname).append((opname, values))
    return opAction

def make_dict_action(argtype=None, postaction='store', first_is_key=True):
    class DictAction(argparse.Action):
        #def __init__(self, option_strings, dest, nargs=None, **kwargs):
        #    if nargs is not None:
        #    raise ValueError("nargs not allowed")
        #    super(DictAction, self).__init__(option_strings, dest, **kwargs)
        def __call__(self, parser, namespace, values, option_string=None):
            if first_is_key: # For the labeled versions
                key = values[0]
                vals = values[1:]
            else: # For the default values
                key = None
                vals = values
            if argtype is not None:
                if (first_is_key and self.nargs == 2) or (not first_is_key and self.nargs == 1):
                    vals = argtype(vals[0]) # avoid having lists in this case
                else:
                    vals = [argtype(val) for val in vals]
            ledict = getattr(namespace, self.dest, {})
            if ledict is None:
                ledict = {}
            if postaction=='store':
                ledict[key] = vals
            elif postaction=='append':
                lelist = ledict.get(key, [])
                lelist.append(vals)
                ledict[key] = lelist
            setattr(namespace, self.dest, ledict)
    return DictAction


class ArgumentProcessingError(Exception):
    pass

class AppendTupleAction(argparse.Action):
    ''' A variation on the 'append' builtin action from argparse, but uses tuples for the internal groupings instead of lists ''' 
    def __call__(self, parser, args, values, option_string=None):
        if (not hasattr(args, self.dest)) or getattr(args, self.dest) is None:
            setattr(args, self.dest, list())
        getattr(args, self.dest).append(tuple(values))

def float_range(string):
    """Tries to parse a string either as one individual float value
    or one of the following patterns:
    
    first:last:increment
    first:last|steps
    first:last
    
    (The last one is equivalent to first:last|50.)
    Returns either float or numpy array.
    """
    try:
        res = float(string)
        return res
    except ValueError:
        import re
        steps = None
        match = re.match(r's?([^:]+):([^|]+)\|(.+)', string)
        if match:
            steps = int(match.group(3))
        else:
            match = re.match(r's?([^:]+):([^:]+):(.+)', string)
            if match:
                increment = float(match.group(3))
            else:
                match = re.match(r's?([^:]+):(.+)', string)
                if match:
                    steps = 50
                else: 
                    argparse.ArgumentTypeError('Invalid float/sequence format: "%s"' % string)
        first = float(match.group(1))
        last = float(match.group(2))
        import numpy as np
        if steps is not None:
            return np.linspace(first, last, num=steps)
        else:
            return np.arange(first, last, increment)

def material_spec(string):
    """Tries to parse a string as a material specification, i.e. a 
    real or complex number or one of the string in built-in Lorentz-Drude models.

    Tries to interpret the string as 1) float, 2) complex, 3) Lorentz-Drude key.
    Raises argparse.ArgumentTypeError on failure.
    """
    from .cymaterials import lorentz_drude
    if string in lorentz_drude.keys():
        return string
    else:
        try: lemat = float(string)
        except ValueError:
            try: lemat = complex(string)
            except ValueError as ve:
                raise argpares.ArgumentTypeError("Material specification must be a supported material name %s, or a number" % (str(lorentz_drude.keys()),)) from ve
        return lemat

class ArgParser:
    ''' Common argument parsing engine for QPMS python CLI scripts. '''
    
    def __add_planewave_argparse_group(ap):
        pwgrp = ap.add_argument_group('Incident wave specification', """
        Incident wave direction is given in terms of ISO polar and azimuthal angles θ, φ, 
        which translate into cartesian coordinates as r̂ = (x, y, z) = (sin(θ) cos(φ), sin(θ) sin(φ), cos(θ)).

        Wave polarisation is given in terms of parameters ψ, χ, where ψ is the angle between a polarisation 
        ellipse axis and meridian tangent θ̂, and tg χ determines axes ratio; 
        the electric field in the origin is then

        E⃗ = cos(χ) (cos(ψ) θ̂ + sin(ψ) φ̂) + i sin(χ) (sin(ψ) θ̂ + cos(ψ) φ̂).

        All the angles are given as multiples of π/2.
        """ # TODO EXAMPLES
        )
        pwgrp.add_argument("-φ", "--phi", type=float, default=0,
            help='Incident wave asimuth in multiples of π/2.')
        pwgrp.add_argument("-θ", "--theta", type=float_range, default=0,
            help='Incident wave polar angle in multiples of π/2. This might be a sequence in format FIRST:LAST:INCREMENT.')
        pwgrp.add_argument("-ψ", "--psi", type=float, default=0,
            help='Angle between polarisation ellipse axis and meridian tangent θ̂ in multiples of π/2.')
        pwgrp.add_argument("-χ", "--chi", type=float, default=0,
            help='Polarisation parameter χ in multiples of π/2. 0 for linear, 0.5 for circular pol.')

    def __add_manyparticle_argparse_group(ap):
        mpgrp = ap.add_argument_group('Many particle specification', "TODO DOC")
        mpgrp.add_argument("-p", "--position", nargs='+', action=make_dict_action(argtype=float, postaction='append',
            first_is_key=False), help="Particle positions, cartesion coordinates (default particle properties)")
        mpgrp.add_argument("+p", "++position", nargs='+', action=make_dict_action(argtype=float, postaction='append', 
            first_is_key=True), help="Particle positions, cartesian coordinates (labeled)")
        mpgrp.add_argument("-L", "--lMax", nargs=1, default={},
                action=make_dict_action(argtype=int, postaction='store', first_is_key=False,),
                help="Cutoff multipole degree (default)")
        mpgrp.add_argument("+L", "++lMax", nargs=2, 
                action=make_dict_action(argtype=int, postaction='store', first_is_key=True,),
            help="Cutoff multipole degree (labeled)")
        mpgrp.add_argument("-m", "--material", nargs=1, default={},
                action=make_dict_action(argtype=material_spec, postaction='store', first_is_key=False,),
                help='particle material (Au, Ag, ... for Lorentz-Drude or number for constant refractive index) (default)')
        mpgrp.add_argument("+m", "++material", nargs=2,
                action=make_dict_action(argtype=material_spec, postaction='store', first_is_key=True,),
                help='particle material (Au, Ag, ... for Lorentz-Drude or number for constant refractive index) (labeled)')
        mpgrp.add_argument("-r", "--radius", nargs=1, default={},
                action=make_dict_action(argtype=float, postaction='store', first_is_key=False,),
                help='particle radius (sphere or cylinder; default)')
        mpgrp.add_argument("+r", "++radius", nargs=2,
                action=make_dict_action(argtype=float, postaction='store', first_is_key=True,),
                help='particle radius (sphere or cylinder; labeled)')
        mpgrp.add_argument("-H", "--height", nargs=1, default={},
                action=make_dict_action(argtype=float, postaction='store', first_is_key=False,),
                help='particle radius (cylinder; default)')
        mpgrp.add_argument("+H", "++height", nargs=2,
                action=make_dict_action(argtype=float, postaction='store', first_is_key=True,),
                help='particle radius (cylinder; labeled)')

    atomic_arguments = {
            'rectlattice2d_periods': lambda ap: ap.add_argument("-p", "--period", type=float, nargs='+', required=True, help='square/rectangular lattice periods', metavar=('px','[py]')),
            'rectlattice2d_counts': lambda ap: ap.add_argument("--size", type=int, nargs=2, required=True, help='rectangular array size (particle column, row count)', metavar=('NCOLS', 'NROWS')),
            'single_frequency_eV': lambda ap: ap.add_argument("-f", "--eV", type=float, required=True, help='radiation angular frequency in eV'),
            'multiple_frequency_eV_optional': lambda ap: ap.add_argument("-f", "--eV", type=float, nargs='*', help='radiation angular frequency in eV (additional)'),
            'seq_frequency_eV': lambda ap: ap.add_argument("-F", "--eV-seq", type=float, nargs=3, required=True, help='uniform radiation angular frequency sequence in eV', metavar=('FIRST', 'INCREMENT', 'LAST')),
            'real_frequencies_eV_ng': lambda ap: ap.add_argument("-f", "--eV", type=float_range, nargs=1, action='append', required=True, help='Angular frequency (or angular frequency range) in eV'), # nargs='+', action='extend' would be better, but action='extend' requires python>=3.8
            'single_material': lambda ap: ap.add_argument("-m", "--material", help='particle material (Au, Ag, ... for Lorentz-Drude or number for constant refractive index)', type=material_spec, required=True),
            'single_radius': lambda ap: ap.add_argument("-r", "--radius", type=float, required=True, help='particle radius (sphere or cylinder)'),
            'single_height': lambda ap: ap.add_argument("-H", "--height", type=float, help='cylindrical particle height; if not provided, particle is assumed to be spherical'),
            'single_kvec2': lambda ap: ap.add_argument("-k", '--kx-lim', nargs=2, type=float, required=True, help='k vector', metavar=('KX_MIN', 'KX_MAX')),
            'kpi': lambda ap: ap.add_argument("--kpi", action='store_true', help="Indicates that the k vector is given in natural units instead of SI, i.e. the arguments given by -k shall be automatically multiplied by pi / period (given by -p argument)"),
            'bg_real_refractive_index': lambda ap: ap.add_argument("-n", "--refractive-index", type=float, default=1., help='background medium strictly real refractive index'),
            'bg_analytical': lambda ap: ap.add_argument("-B", "--background", type=material_spec, default=1., help="Background medium specification (constant real or complex refractive index, or supported material label)"),
            'single_lMax': lambda ap: ap.add_argument("-L", "--lMax", type=int, required=True, default=3, help='multipole degree cutoff'),
            'single_lMax_extend': lambda ap: ap.add_argument("--lMax-extend", type=int, required=False, default=6, help='multipole degree cutoff for T-matrix calculation (cylindrical particles only'),
            'outfile': lambda ap: ap.add_argument("-o", "--output", type=str, required=False, help='output path (if not provided, will be generated automatically)'), # TODO consider type=argparse.FileType('w')
            'plot_out': lambda ap: ap.add_argument("-O", "--plot-out", type=str, required=False, help="path to plot output (optional)"),
            'plot_do': lambda ap: ap.add_argument("-P", "--plot", action='store_true', help="if -p not given, plot to a default path"),
            'lattice2d_basis': lambda ap: ap.add_argument("-b", "--basis-vector", nargs='+', action=AppendTupleAction, help="basis vector in xy-cartesian coordinates (two required)", required=True, dest='basis_vectors', metavar=('X', 'Y')),
            'planewave_pol_angles': __add_planewave_argparse_group,
            'multi_particle': __add_manyparticle_argparse_group,
    }

    feature_sets_available = { # name : (description, dependencies, atoms not in other dependencies, methods called after parsing, "virtual" features provided) 
            'const_real_background': ("Background medium with constant real refractive index", (), ('bg_real_refractive_index',), ('_eval_const_background_epsmu',), ('background', 'background_analytical')),
            'background' : ("Most general background medium specification currently supported", ('background_analytical',), (), (), ()),
            'background_analytical' : ("Background medium model holomorphic for 'reasonably large' complex frequency areas", (), ('bg_analytical',), ('_eval_analytical_background_epsmugen',), ('background',)),
            'single_particle': ("Single particle definition (shape [currently spherical or cylindrical]) and materials, incl. background)", ('background',), ('single_material', 'single_radius', 'single_height', 'single_lMax_extend'), ('_eval_single_tmgen',), ()),
            'multi_particle': ("One or more particle definition (shape [curently spherical or cylindrical]), materials, and positions)", ('background',), ('multi_particle',), ('_process_multi_particle',), ()),
            'single_lMax': ("Single particle lMax definition", (), ('single_lMax',), (), ()),
            'single_omega': ("Single angular frequency", (), ('single_frequency_eV',), ('_eval_single_omega',), ()),
            'omega_seq': ("Equidistant real frequency range with possibility of adding individual frequencies", (), ('seq_frequency_eV', 'multiple_frequency_eV_optional',), ('_eval_omega_seq',), ()),
            'omega_seq_real_ng': ("Equidistant real frequency ranges or individual frequencies (new syntax)", (), ('real_frequencies_eV_ng',), ('_eval_omega_seq_real_ng',), ()),
            'lattice2d': ("Specification of a generic 2d lattice (spanned by the x,y axes)", (), ('lattice2d_basis',), ('_eval_lattice2d',), ()),
            'rectlattice2d': ("Specification of a rectangular 2d lattice; conflicts with lattice2d", (), ('rectlattice2d_periods',), ('_eval_rectlattice2d',), ()),
            'rectlattice2d_finite': ("Specification of a rectangular 2d lattice; conflicts with lattice2d", ('rectlattice2d',), ('rectlattice2d_counts',), (), ()),
            'planewave': ("Specification of a normalised plane wave (typically used for scattering) with a full polarisation state", (), ('planewave_pol_angles',), ("_process_planewave_angles",), ()),
    }


    def __init__(self, features=[]):
        prefix_chars = '+-' if 'multi_particle' in features else '-'
        self.ap = argparse.ArgumentParser(prefix_chars=prefix_chars)
        self.features_enabled = set()
        self.call_at_parse_list = []
        self.parsed = False
        for feat in features:
            self.add_feature(feat)
        self._emg_register = {} # EpsMuGenerator dictionary to avoid recreating equivalent instances; filled by _add_emg()
        self._tmg_register = {} # TMatrixGenerator dictionary to avoid recreating equivalent instances; filled by _add_tmg()
        self._bspec_register = {} # Dictionary of used BaseSpecs to keep the equivalent instances unique; filled by _add_bspec()

    def _add_emg(self, emgspec):
        """Looks up whether if an EpsMuGenerator from given material_spec has been already registered, and if not, creates a new one"""
        from .cymaterials import EpsMu, EpsMuGenerator, lorentz_drude
        if emgspec in self._emg_register.keys():
            return self._emg_register[emgspec]
        else:
            if isinstance(emgspec, (float, complex)):
                emg = EpsMuGenerator(EpsMu(emgspec**2))
            else:
                emg = EpsMuGenerator(lorentz_drude[emgspec])
            self._emg_register[emgspec] = emg
            return emg

    def _add_tmg(self, tmgspec):
        """Looks up whether if a T-matrix from given T-matrix specification tuple has been already registered, and if not, creates a new one
        
        T-matrix specification shall be of the form 
        (bg_material_spec, fg_material_spec, shape_spec) where shape_spec is
        (radius, height, lMax_extend)
        """
        if tmgspec in self._tmg_register.keys():
            return self._tmg_register[tmgspec]
        else:
            from .cytmatrices import TMatrixGenerator
            bgspec, fgspec, (radius, height, lMax_extend) = tmgspec
            bg = self._add_emg(bgspec)
            fg = self._add_emg(fgspec)
            if height is None:
                tmgen = TMatrixGenerator.sphere(bg, fg, radius)
            else:
                tmgen = TMatrixGenerator.cylinder(bg, fg, radius, height, lMax_extend=lMax_extend)
            self._tmg_register[tmgspec] = tmgen
            return tmgen

    def _add_bspec(self, key):
        if key in self._bspec_register.keys():
            return self._bspec_register[key]
        else:
            from .cybspec import BaseSpec
            if isinstance(key, BaseSpec):
                bspec = key
            elif isinstance(key, int):
                bspec = self._add_bspec(BaseSpec(lMax=key))
            else: raise TypeError("Can't register this as a BaseSpec")
            self._bspec_register[key] = bspec
            return bspec

    def add_feature(self, feat):
        if feat not in self.features_enabled:
            if feat not in ArgParser.feature_sets_available:
                raise ValueError("Unknown ArgParser feature: %s" % feat)
            #resolve dependencies
            _, deps, atoms, atparse, provides_virtual = ArgParser.feature_sets_available[feat]
            for dep in deps:
                self.add_feature(dep)
            for atom in atoms: # maybe check whether that atom has already been added sometimes in the future?
                ArgParser.atomic_arguments[atom](self.ap)
            for methodname in atparse:
                self.call_at_parse_list.append(methodname)
            self.features_enabled.add(feat)
            for feat_virt in provides_virtual:
                self.features_enabled.add(feat_virt)

    def add_argument(self, *args, **kwargs):
        '''Add a custom argument directly to the standard library ArgParser object'''
        self.ap.add_argument(*args, **kwargs)

    def parse_args(self, process_data = True,  *args, **kwargs):
        self.args = self.ap.parse_args(*args, **kwargs)
        if process_data:
            for method in self.call_at_parse_list:
                try:
                    getattr(self, method)()
                except ArgumentProcessingError:
                    err = sys.exc_info()[1]
                    self.ap.error(str(err))
        return self.args
    
    def __getattr__(self, name):
        return getattr(self.args, name)
        

    # Methods to initialise the related data structures:

    def _eval_const_background_epsmu(self): # feature: const_real_background
        self.args.background = self.args.refractive_index
        self._eval_analytical_background_epsmugen()
    
    def _eval_analytical_background_epsmugen(self): # feature: background_analytical
        a = self.args
        from .cymaterials import EpsMu
        if isinstance(a.background, (float, complex)):
            self.background_epsmu = EpsMu(a.background**2)
        self.background_emg = self._add_emg(a.background)

    def _eval_single_tmgen(self): # feature: single_particle
        a = self.args
        from .cymaterials import EpsMuGenerator, lorentz_drude
        from .cytmatrices import TMatrixGenerator
        self.foreground_emg = self._add_emg(a.material)
        self.tmgen = self._add_tmg((a.background, a.material, (a.radius, a.height, a.lMax_extend)))
        self.bspec = self._add_bspec(a.lMax)

    def _eval_single_omega(self): # feature: single_omega
        from .constants import eV, hbar
        self.omega = self.args.eV * eV / hbar

    def _eval_omega_seq(self): # feature: omega_seq
        import numpy as np
        from .constants import eV, hbar
        start, step, stop = self.args.eV_seq
        self.omegas = np.arange(start, stop, step)
        if self.args.eV:
            self.omegas = np.concatenate((self.omegas, np.array(self.args.eV)))
            self.omegas.sort()
        self.omegas *= eV/hbar

    def _eval_omega_seq_real_ng(self): # feature: omega_seq_real_ng
        import numpy as np
        from .constants import eV, hbar
        eh = eV / hbar
        self.omegas = [omega_eV * eh for  omega_eV in flatten(self.args.eV)]
        self.omega_max = max(om if isinstance(om, float) else max(om) for om in self.omegas)
        self.omega_min = min(om if isinstance(om, float) else min(om) for om in self.omegas)
        self.omega_singles = [om for om in self.omegas if isinstance(om, float)]
        self.omega_ranges = [om for om in self.omegas if not isinstance(om, float)]
        self.omega_descr = ("%geV" % (self.omega_max / eh)) if (self.omega_max == self.omega_min) else (
                "%g–%geV" % (self.omega_min / eh, self.omega_max / eh))
        self.allomegas = []
        for om in self.omegas:
            if isinstance(om, float):
                self.allomegas.append(om)
            else:
                self.allomegas.extend(om)
        self.allomegas = np.unique(self.allomegas)


    def _eval_lattice2d(self): # feature: lattice2d
        l = len(self.args.basis_vectors)
        if l != 2: raise ValueError('Two basis vectors must be specified (have %d)' % l)
        from .qpms_c import lll_reduce
        self.direct_basis = lll_reduce(self.args.basis_vectors, delta=1.)
        import numpy as np
        self.reciprocal_basis1 = np.linalg.inv(self.direct_basis)
        self.reciprocal_basis2pi = 2 * np.pi * self.reciprocal_basis1
    
    def _eval_rectlattice2d(self): # feature: rectlattice2d
        a = self.args
        l = len(a.period)
        if (l == 1): # square lattice
            a.period = (a.period[0], a.period[0])
        else:
            a.period = (a.period[0], a.period[1])
            if (l > 2):
                raise ValueError("At most two lattice periods allowed for a rectangular lattice (got %d)" % l)

        import numpy as np
        a.basis_vectors = [(a.period[0], 0.), (0., a.period[1])]
        self.direct_basis = np.array(a.basis_vectors)
        self.reciprocal_basis1 = np.linalg.inv(self.direct_basis)
        self.reciprocal_basis2pi = 2 * np.pi * self.reciprocal_basis1

    def _process_planewave_angles(self): #feature: planewave
        import math
        pi2 = math.pi/2
        a = self.args
        a.chi = a.chi * pi2
        a.psi = a.psi * pi2
        a.theta = a.theta * pi2
        a.phi = a.phi * pi2

    def _process_multi_particle(self): # feature: multi_particle
        a = self.args
        self.tmspecs = {}
        self.tmgens = {}
        self.bspecs = {}
        self.positions = {}
        pos13, pos23, pos33 = False, False, False # used to 
        if len(a.position.keys()) == 0:
            warnings.warn("No particle position (-p or +p) specified, assuming single particle in the origin / single particle per unit cell!")
            a.position[None] = [(0.,0.,0.)]
        for poslabel in a.position.keys():
            try:
                lMax = a.lMax.get(poslabel, False) or a.lMax[None]
                radius = a.radius.get(poslabel, False) or a.radius[None]
                # Height is "inherited" only together with radius
                height = a.height.get(poslabel, None) if poslabel in a.radius.keys() else a.height.get(None, None)
                if hasattr(a, 'lMax_extend'):
                    lMax_extend = a.lMax_extend.get(poslabel, False) or a.lMax_extend.get(None, False) or None
                else:
                    lMax_extend = None
                material = a.material.get(poslabel, False) or a.material[None]
            except (TypeError, KeyError) as exc:
                if poslabel is None:
                    raise ArgumentProcessingError("Unlabeled particles' positions (-p) specified, but some default particle properties are missing (--lMax, --radius, and --material have to be specified)") from exc
                else:
                    raise ArgumentProcessingError(("Incomplete specification of '%s'-labeled particles: you must"
                        "provide at least ++lMax, ++radius, ++material arguments with the label, or the fallback arguments"
                        "--lMax, --radius, --material.")%(str(poslabel),)) from exc
            tmspec = (a.background, material, (radius, height, lMax_extend)) 
            self.tmspecs[poslabel] = tmspec
            self.tmgens[poslabel] = self._add_tmg(tmspec)
            self.bspecs[poslabel] = self._add_bspec(lMax)
            poslist_cured = []
            for pos in a.position[poslabel]:
                if len(pos) == 1:
                    pos_cured = (0., 0., pos[0])
                    pos13 = True
                elif len(pos) == 2:
                    pos_cured = (pos[0], pos[1], 0.)
                    pos23 = True
                elif len(pos) == 3:
                    pos_cured = pos
                    pos33 = True
                else:
                    raise argparse.ArgumentTypeError("Each -p / +p argument requires 1 to 3 cartesian coordinates")
                poslist_cured.append(pos_cured)
            self.positions[poslabel] = poslist_cured
        if pos13 and pos23:
            warnings.warn("Both 1D and 2D position specifications used. The former are interpreted as z coordinates while the latter as x, y coordinates")

    def get_particles(self):
        """Creates a list of Particle instances that can be directly used in ScatteringSystem.create().

        Assumes that self._process_multi_particle() has been already called.
        """
        from .qpms_c import Particle
        plist = []
        for poslabel, poss in self.positions.items():
            t = self.tmgens[poslabel]
            bspec = self.bspecs[poslabel]
            plist.extend([Particle(pos, t, bspec=bspec) for pos in poss])
        return plist

