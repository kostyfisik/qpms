'''
Common snippets for argument processing in command line scripts; legacy scripts use scripts_common.py instead.
'''

import argparse

class ArgParser:
    ''' Common argument parsing engine for QPMS python CLI scripts. '''
    atomic_arguments = {
            'sqlat_period': lambda ap: ap.add_argument("-p", "--period", type=float, required=True, help='square lattice period'),
            'rectlat_Nx': lambda ap: ap.add_argument("--Nx", type=int, required=True, help='array size x'),
            'rectlat_Ny': lambda ap: ap.add_argument("--Ny", type=int, required=True, help='array size y'),
            'single_frequency_eV': lambda ap: ap.add_argument("-f", "--eV", type=float, required=True, help='radiation angular frequency in eV'),
            'single_material': lambda ap: ap.add_argument("-m", "--material", help='particle material (Au, Ag, ... for Lorentz-Drude or number for constant refractive index)', default='Au', required=True),
            'single_radius': lambda ap: ap.add_argument("-r", "--radius", type=float, required=True, help='particle radius (sphere or cylinder)'),
            'single_height': lambda ap: ap.add_argument("-H", "--height", type=float, help='cylindrical particle height; if not provided, particle is assumed to be spherical'),
            'single_kvec2': lambda ap: ap.add_argument("-k", '--kx-lim', nargs=2, type=float, required=True, help='k vector', metavar=('KX_MIN', 'KX_MAX')),
            'kpi': lambda ap: ap.add_argument("--kpi", action='store_true', help="Indicates that the k vector is given in natural units instead of SI, i.e. the arguments given by -k shall be automatically multiplied by pi / period (given by -p argument)"),
            'bg_refractive_index': lambda ap: ap.add_argument("-n", "--refractive-index", type=float, default=1.52, help='background medium refractive index'),
            'single_lMax': lambda ap: ap.add_argument("-L", "--lMax", type=int, required=True, default=3, help='multipole degree cutoff'),
            'single_lMax_extend': lambda ap: ap.add_argument("--lMax-extend", type=int, required=False, default=6, help='multipole degree cutoff for T-matrix calculation (cylindrical particles only'),
            'outfile': lambda ap: ap.add_argument("-o", "--output", type=str, required=False, help='output path (if not provided, will be generated automatically)'),
            'plot_out': lambda ap: ap.add_argument("-O", "--plot-out", type=str, required=False, help="path to plot output (optional)"),
            'plot_do': lambda ap: ap.add_argument("-P", "--plot", action='store_true', help="if -p not given, plot to a default path"),
    }

    feature_sets_available = { # name : (description, dependencies, atoms not in other dependencies, methods called after parsing) 
            'background': ("Background medium definition (currently only constant epsilon supported)", (), ('bg_refractive_index',), ('_eval_background_epsmu',)),
            'single_particle': ("Single particle definition (shape [currently spherical or cylindrical]) and materials, incl. background)", ('background',), ('single_material', 'single_radius', 'single_height', 'single_lMax_extend'), ('_eval_single_tmgen',)),
            'single_lMax': ("Single particle lMax definition", (), ('single_lMax',), ()),
            'single_omega': ("Single angular frequency", (), ('single_frequency_eV',), ('_eval_single_omega',)),
    }


    def __init__(self, features=[]):
        self.ap = argparse.ArgumentParser()
        self.features_enabled = set()
        self.call_at_parse_list = []
        self.parsed = False
        for feat in features:
            self.add_feature(feat)

    def add_feature(self, feat):
        if feat not in self.features_enabled:
            if feat not in ArgParser.feature_sets_available:
                raise ValueError("Unknown ArgParser feature: %s", feat)
            #resolve dependencies
            _, deps, atoms, atparse = ArgParser.feature_sets_available[feat]
            for dep in deps:
                self.add_feature(dep)
            for atom in atoms: # maybe check whether that atom has already been added sometimes in the future?
                ArgParser.atomic_arguments[atom](self.ap)
            for methodname in atparse:
                self.call_at_parse_list.append(methodname)
            self.features_enabled.add(feat)

    def add_argument(self, *args, **kwargs):
        '''Add a custom argument directly to the standard library ArgParser object'''
        self.ap.add_argument(*args, **kwargs)

    def parse_args(self, process_data = True,  *args, **kwargs):
        self.args = self.ap.parse_args(*args, **kwargs)
        if process_data:
            for method in self.call_at_parse_list:
                getattr(self, method)()
        return self.args
    
    def __getattr__(self, name):
        return getattr(self.args, name)
        

    # Methods to initialise the related data structures:

    def _eval_background_epsmu(self): # feature: background
        from .cymaterials import EpsMu
        self.background_epsmu = EpsMu(self.args.refractive_index**2)

    def _eval_single_tmgen(self): # feature: single_particle
        a = self.args
        from .cymaterials import EpsMuGenerator, lorentz_drude
        from .cytmatrices import TMatrixGenerator
        if a.material in lorentz_drude.keys():
            self.foreground_emg = EpsMuGenerator(lorentz_drude[a.material])
        else:
            try: lemat = float(a.material)
            except ValueError:
                try: lemat = complex(a.material)
                except ValueError as ve:
                    raise ValueError("--material must be either a label such as 'Ag', 'Au', or a number") from ve
            a.material = lemat
            self.foreground_emg = EpsMuGenerator(EpsMu(a.material**2))

        if a.height is None:
            self.tmgen = TMatrixGenerator.sphere(self.background_epsmu, self.foreground_emg, a.radius)
        else:
            self.tmgen = TMatrixGenerator.cylinder(self.background_epsmu, self.foreground_emg, a.radius, a.height, lMax_extend = a.lMax_extend)

    def _eval_single_omega(self): # feature: single_omega
        from .constants import eV, hbar
        self.omega = self.args.eV * eV / hbar


