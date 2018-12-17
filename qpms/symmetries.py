from sympy.combinatorics import Permutation, PermutationGroup
Permutation.print_cyclic = True
import cmath
from cmath import exp, pi
from math import sqrt
import numpy as np
np.set_printoptions(linewidth=200)
import qpms
import numbers
ň = None

def grouprep_try(tdict, src, im, srcgens, imgens, immultop = None, imcmp = None):
    tdict[src] = im
    for i in range(len(srcgens)):
        new_src = src * srcgens[i]
        new_im = (im * imgens[i]) if (immultop is None) else immultop(im, imgens[i])
        if new_src not in tdict.keys():
            grouprep_try(tdict, new_src, new_im, srcgens, imgens, immultop, imcmp)
        elif ((new_im != tdict[new_src]) if (imcmp is None) else (not imcmp(new_im, tdict[new_src]))): # check consistency
            print(src, ' * ', srcgens[i], ' --> ', new_src)
            print(im)
            print(' * ')
            print(imgens[i])
            print(' --> ')
            print(new_im)
            print(' != ')
            print(tdict[new_src])
            raise ValueError("Homomorphism inconsistency detected")                
    return

class SVWFPointGroupInfo: # only for point groups, coz in svwf_rep() I use I_tyty, not I_ptypty or something alike
    def __init__(self,
                 name,
                 permgroupgens, # permutation group generators
                 irrepgens_dict, # dictionary with irrep generators,
                 svwf_rep_gen_func, # function that generates a tuple with svwf representation generators
                ):
        self.name = name
        self.permgroupgens = permgroupgens
        self.permgroup = PermutationGroup(*permgroupgens)
        self.irrepgens_dict = irrepgens_dict
        self.svwf_rep_gen_func = svwf_rep_gen_func
        self.irreps = dict()
        for irrepname, irrepgens in irrepgens_dict.items():
            is1d = isinstance(irrepgens[0], int)
            irrepdim = 1 if is1d else irrepgens[0].shape[0]
            self.irreps[irrepname] = generate_grouprep(self.permgroup, 
                                                       1 if is1d else np.eye(irrepdim),
                                                       permgroupgens, irrepgens,
                                                       immultop = None if is1d else np.dot,
                                                       imcmp = None if is1d else np.allclose
                                                      )
    
    def svwf_rep(self, lMax, *rep_gen_func_args, **rep_gen_func_kwargs):
        '''
        This method generates full SVWF (reducible) representation of the group.
        '''
        svwfgens = self.svwf_rep_gen_func(lMax, *rep_gen_func_args, **rep_gen_func_kwargs)
        my, ny = qpms.get_mn_y(lMax)
        nelem = len(my)
        I_tyty = np.moveaxis(np.eye(2)[:,:,ň,ň] * np.eye(nelem), 2,1)
        return generate_grouprep(self.permgroup, I_tyty, self.permgroupgens, svwfgens, immultop = mmult_tyty, imcmp = np.allclose)
    
    def svwf_irrep_projectors(self, lMax, *rep_gen_func_args, **rep_gen_func_kwargs):
        return gen_point_group_svwfrep_projectors(self.permgroup, self.irreps, self.svwf_rep(lMax, *rep_gen_func_args, **rep_gen_func_kwargs))

# srcgroup is expected to be PermutationGroup and srcgens of the TODO
# imcmp returns True if two elements of the image group are 'equal', otherwise False
def generate_grouprep(srcgroup, im_identity, srcgens, imgens, immultop = None, imcmp = None):
    sz = srcgens[0].size
    for g in srcgens:
        if g.size != sz:
            raise ValueError('All the generators must have the same "size"')
    tdict = dict()
    grouprep_try(tdict, Permutation(sz-1), im_identity, srcgens, imgens, immultop = immultop, imcmp = imcmp)
    if(srcgroup.order() != len(tdict.keys())): # basic check
        raise ValueError('The supplied "generators" failed to generate the preimage group: ', 
                         srcgroup.order(), " != ", len(tdict.keys()))
    return tdict
    
# matrices appearing in 2d representations of common groups as used in Bradley, Cracknell p. 61 (with arabic names instead of greek, because lambda is a keyword)
epsilon = np.eye(2)
alif = np.array(((-1/2,-sqrt(3)/2),(sqrt(3)/2,-1/2)))
bih = np.array(((-1/2,sqrt(3)/2),(-sqrt(3)/2,-1/2)))
kaf = np.array(((0,1),(1,0)))
lam = np.array(((1,0),(0,-1)))
ra = np.array(((0,-1),(1,0)))
mim =  np.array(((-1/2,-sqrt(3)/2),(-sqrt(3)/2,1/2)))
nun =  np.array(((-1/2,sqrt(3)/2),(sqrt(3)/2,1/2)))




def mmult_tyty(a, b):
        return(qpms.apply_ndmatrix_left(a, b, (-4,-3)))
def mmult_ptypty(a, b):
    return(qpms.apply_ndmatrix_left(a, b, (-6,-5,-4)))

def gen_point_group_svwfrep_projectors(permgroup, matrix_irreps_dict, sphrep_full):
    '''
    Gives the projection operators $P_kl('\Gamma')$ from Dresselhaus (4.28)
    for all irreps $\Gamma$ of D3h.;
    as an array with indices [k,l,t,y,t,y]
    
    Example of creating last argument:
    sphrep_full = generate_grouprep(D3h_permgroup, I_tyty, D3h_srcgens, [C3_tyty, vfl_tyty, zfl_tyty], 
                           immultop = mmult_tyty, imcmp = np.allclose)
    '''
    order = permgroup.order()
    sphreps = dict()
    nelem = sphrep_full[permgroup[0]].shape[-1] # quite ugly hack
    for repkey, matrixrep in matrix_irreps_dict.items():
        arepmatrix = matrixrep[permgroup[0]] # just one of the matrices to get the shape etc
        if isinstance(arepmatrix, numbers.Number):
            dim = 1 # repre dimension
            preprocess = lambda x: np.array([[x]])
        elif isinstance(arepmatrix, np.ndarray):
            if(len(arepmatrix.shape)) != 2 or arepmatrix.shape[0] != arepmatrix.shape[1]:
                raise ValueError("Arrays representing irrep matrices must be of square shape")
            dim = arepmatrix.shape[0]
            preprocess = lambda x: x
        else: 
            raise ValueError("Irrep is not a square array or number")        
        sphrep = np.zeros((dim,dim,2,nelem,2,nelem), dtype=complex)
        for i in permgroup.elements:
            sphrep += preprocess(matrixrep[i]).conj().transpose()[:,:,ň,ň,ň,ň] * sphrep_full[i]
        sphrep *= dim / order
        # clean the nonexact values here 
        for x in [0, 0.5, -0.5, 0.5j, -0.5j]:
            sphrep[np.isclose(sphrep,x)]=x
        sphreps[repkey] = sphrep
    return sphreps

# Group D3h; mostly legacy code (kept because of the the honeycomb lattice K-point code, whose generalised version not yet implemented)
# Note that the size argument of permutations is necessary, otherwise e.g. c*c and  b*b would not be evaluated equal
# N.B. the weird elements as Permutation(N) – it means identity permutation of size N+1.
rot3_perm = Permutation(0,1,2, size=5) # C3 rotation
xflip_perm = Permutation(0,2, size=5) # vertical mirror
zflip_perm = Permutation(3,4, size=5) # horizontal mirror
D3h_srcgens = [rot3_perm,xflip_perm,zflip_perm]
D3h_permgroup = PermutationGroup(*D3h_srcgens) # D3h

D3h_irreps = {
    # Bradley, Cracknell p. 61
    "E'" : generate_grouprep(D3h_permgroup, epsilon, D3h_srcgens, [alif, lam, epsilon], immultop = np.dot, imcmp = np.allclose),
    "E''" : generate_grouprep(D3h_permgroup, epsilon, D3h_srcgens, [alif, lam, -epsilon], immultop = np.dot, imcmp = np.allclose),
    # Bradley, Cracknell p. 59, or Dresselhaus, Table A.14 (p. 482)
    "A1'" : generate_grouprep(D3h_permgroup, 1, D3h_srcgens, [1,1,1]),
    "A2'" : generate_grouprep(D3h_permgroup, 1, D3h_srcgens, [1,-1,1]),
    "A1''" : generate_grouprep(D3h_permgroup, 1, D3h_srcgens, [1,-1,-1]),
    "A2''" : generate_grouprep(D3h_permgroup, 1, D3h_srcgens, [1,1,-1]),
}

#TODO lepší název fce; legacy, use group_info['D3h'].generate_grouprep() instead
def gen_point_D3h_svwf_rep(lMax, vflip = 'x'):
    '''
    Gives the projection operators $P_kl('\Gamma')$ from Dresselhaus (4.28)
    for all irreps $\Gamma$ of D3h.;
    as an array with indices [k,l,t,y,t,y]
    '''

    my, ny = qpms.get_mn_y(lMax)
    nelem = len(my)
    C3_yy = qpms.WignerD_yy_fromvector(lMax, np.array([0,0,2*pi/3]))
    C3_tyty = np.moveaxis(np.eye(2)[:,:,ň,ň] * C3_yy, 2,1)
    zfl_tyty = qpms.zflip_tyty(lMax)
    #yfl_tyty = qpms.yflip_tyty(lMax)
    #xfl_tyty = qpms.xflip_tyty(lMax)
    vfl_tyty = qpms.yflip_tyty(lMax) if vflip == 'y' else qpms.xflip_tyty(lMax)
    I_tyty = np.moveaxis(np.eye(2)[:,:,ň,ň] * np.eye(nelem), 2,1)
    order = D3h_permgroup.order()
    sphrep_full = generate_grouprep(D3h_permgroup, I_tyty, D3h_srcgens, [C3_tyty, vfl_tyty, zfl_tyty], 
                           immultop = mmult_tyty, imcmp = np.allclose)
    sphreps = dict()
    for repkey, matrixrep in D3h_irreps.items():
        arepmatrix = matrixrep[rot3_perm] # just one of the matrices to get the shape etc
        if isinstance(arepmatrix, numbers.Number):
            dim = 1 # repre dimension
            preprocess = lambda x: np.array([[x]])
        elif isinstance(arepmatrix, np.ndarray):
            if(len(arepmatrix.shape)) != 2 or arepmatrix.shape[0] != arepmatrix.shape[1]:
                raise ValueError("Arrays representing irrep matrices must be of square shape")
            dim = arepmatrix.shape[0]
            preprocess = lambda x: x
        else: 
            raise ValueError("Irrep is not a square array or number")
        sphrep = np.zeros((dim,dim,2,nelem,2,nelem), dtype=complex)
        for i in D3h_permgroup.elements:
            sphrep += preprocess(matrixrep[i]).conj().transpose()[:,:,ň,ň,ň,ň] * sphrep_full[i]
        sphrep *= dim / order
        # clean the nonexact values here 
        for x in [0, 0.5, -0.5, 0.5j, -0.5j]:
            sphrep[np.isclose(sphrep,x)]=x
        sphreps[repkey] = sphrep
    return sphreps
        
def gen_hexlattice_Kpoint_svwf_rep(lMax, psi, vflip = 'x'):
    my, ny = qpms.get_mn_y(lMax)
    nelem = len(my)
    C3_yy = qpms.WignerD_yy_fromvector(lMax, np.array([0,0,2*pi/3]))
    C3_tyty = np.moveaxis(np.eye(2)[:,:,ň,ň] * C3_yy, 2,1)
    zfl_tyty = qpms.zflip_tyty(lMax)
    #yfl_tyty = qpms.yflip_tyty(lMax)
    #xfl_tyty = qpms.xflip_tyty(lMax)
    vfl_tyty = qpms.yflip_tyty(lMax) if vflip == 'y' else qpms.xflip_tyty(lMax)
    I_tyty = np.moveaxis(np.eye(2)[:,:,ň,ň] * np.eye(nelem), 2,1)
    hex_C3_K_ptypty = np.diag([exp(-psi*1j*2*pi/3),exp(+psi*1j*2*pi/3)])[:,ň,ň,:,ň,ň] * C3_tyty[ň,:,:,ň,:,:]
    hex_zfl_ptypty = np.eye(2)[:,ň,ň,:,ň,ň] * zfl_tyty[ň,:,:,ň,:,:]
    #hex_xfl_ptypty = np.array([[0,1],[1,0]])[:,ň,ň,:,ň,ň] * xfl_tyty[ň,:,:,ň,:,:]
    hex_vfl_ptypty = np.array([[0,1],[1,0]])[:,ň,ň,:,ň,ň] * vfl_tyty[ň,:,:,ň,:,:]
    hex_I_ptypty = np.eye((2*2*nelem)).reshape((2,2,nelem,2,2,nelem))
    order = D3h_permgroup.order()
    hex_K_sphrep_full = generate_grouprep(D3h_permgroup, hex_I_ptypty, D3h_srcgens, [hex_C3_K_ptypty, hex_vfl_ptypty, hex_zfl_ptypty], 
                           immultop = mmult_ptypty, imcmp = np.allclose)
    hex_K_sphreps = dict()
    for repkey, matrixrep in D3h_irreps.items():
        arepmatrix = matrixrep[rot3_perm] # just one of the matrices to get the shape etc
        if isinstance(arepmatrix, numbers.Number):
            dim = 1 # repre dimension
            preprocess = lambda x: np.array([[x]])
        elif isinstance(arepmatrix, np.ndarray):
            if(len(arepmatrix.shape)) != 2 or arepmatrix.shape[0] != arepmatrix.shape[1]:
                raise ValueError("Arrays representing irrep matrices must be of square shape")
            dim = arepmatrix.shape[0]
            preprocess = lambda x: x
        else: 
            raise ValueError("Irrep is not a square array or number")
        sphrep = np.zeros((dim,dim,2,2,nelem,2,2,nelem), dtype=complex)
        for i in D3h_permgroup.elements:
            sphrep += preprocess(matrixrep[i]).conj().transpose()[:,:,ň,ň,ň,ň,ň,ň] * hex_K_sphrep_full[i]
        sphrep *= dim / order
        # clean the nonexact values here 
        for x in [0, 0.5, -0.5, 0.5j, -0.5j]:
            sphrep[np.isclose(sphrep,x)]=x
        hex_K_sphreps[repkey] = sphrep
    return hex_K_sphreps        

def normalize(v):
    norm = np.linalg.norm(v.reshape((np.prod(v.shape),)), ord=2)
    if norm == 0: 
       return v*np.nan
    return v / norm

def gen_hexlattice_Kpoint_svwf_rep_projectors(lMax, psi, vflip='x', do_bases=False):
    nelem = lMax * (lMax+2)
    projectors = dict()
    if do_bases:
        bases = dict()
    for repi, W in gen_hexlattice_Kpoint_svwf_rep(lMax,psi,vflip=vflip).items():
        totalvecs = 0
        tmplist = list()
        for p in (0,1):
         for t in (0,1):
          for y in range(nelem):
           for ai in range(W.shape[0]):
            for bi in range(W.shape[1]):
                v = np.zeros((2,2,nelem))
                v[p,t,y] = 1
                #v = np.ones((2,2,nelem))
                v1 = np.tensordot(W[ai,bi],v, axes = ([-3,-2,-1],[0,1,2])) 


                if not np.allclose(v1,0):
                    v1 = normalize(v1)
                    for v2 in tmplist:
                        dot = np.tensordot(v1.conjugate(),v2,axes = ([-3,-2,-1],[0,1,2]))
                        if not np.allclose(dot,0):
                            if not np.allclose(np.abs(dot),1):
                                raise ValueError('You have to fix this piece of code.')# TODO maybe I should make sure that the absolute value is around 1
                            break
                    else:
                        totalvecs += 1
                        tmplist.append(v1)
                        #for index, x in np.ndenumerate(v1):
                        #    if x!=0:
                        #        print(index, x)
                        #print('----------')
        theprojector = np.zeros((2,2,nelem,2,2,nelem), dtype = float) 
        if do_bases:
            thebasis = np.zeros((len(tmplist), 2,2,nelem), dtype=complex)
            for i, v in enumerate(tmplist):
                thebasis[i] = v
            bases[repi] = thebasis
        for v in tmplist:
            theprojector += (v[:,:,:,ň,ň,ň] * v.conjugate()[ň,ň,ň,:,:,:]).real # TODO check is it possible to have imaginary elements?
        for x in [0, 1, -1,sqrt(0.5),-sqrt(0.5),0.5,-0.5]:
            theprojector[np.isclose(theprojector,x)]=x
        projectors[repi] = theprojector
    if do_bases:
        return projectors, bases
    else:
        return projectors



point_group_info = { # representation info of some useful point groups
    'C2v' : SVWFPointGroupInfo('C2v',
                 # permutation group generators
                               (Permutation(0,1, size=4)(2,3), # x -> - x mirror operation (i.e. yz mirror plane)
                                Permutation(0,3, size=4)(1,2), # y -> - y mirror operation (i.e. xz mirror plane)
                               ), 
                 # dictionary with irrep generators
                               {
                                    # Bradley, Cracknell p. 58; not sure about the labels / axes here
                                    'A1': (1,1),
                                    'B2': (-1,1),
                                    'A2': (-1,-1),
                                    'B1': (1,-1),    
                               },
                 # function that generates a tuple with svwf representation generators
                               lambda lMax : (qpms.xflip_tyty(lMax), qpms.yflip_tyty(lMax))
    ),
    'D2h' : SVWFPointGroupInfo('D2h',
                 # permutation group generators
                               (Permutation(0,1, size=6)(2,3), # x -> - x mirror operation (i.e. yz mirror plane)
                                Permutation(0,3, size=6)(1,2), # y -> - y mirror operation (i.e. xz mirror plane)
                                Permutation(4,5, size=6)       # z -> - z mirror operation (i.e. xy mirror plane)
                               ), 
                 # dictionary with irrep generators
                               {
                                    # Product of C2v and zflip; not sure about the labels / axes here
                                    "A1'": (1,1,1),
                                    "B2'": (-1,1,1),
                                    "A2'": (-1,-1,1),
                                    "B1'": (1,-1,1),
                                    "A1''": (-1,-1,-1),
                                    "B2''": (1,-1,-1),
                                    "A2''": (1,1,-1),
                                    "B1''": (-1,1,-1),
                               },
                 # function that generates a tuple with svwf representation generators
                               lambda lMax : (qpms.xflip_tyty(lMax), qpms.yflip_tyty(lMax), qpms.zflip_tyty(lMax))
    ),
    'C4v' : SVWFPointGroupInfo('C4v',
                 # permutation group generators
                               (Permutation(0,1,2,3, size=4), #C4 rotation
                                Permutation(0,1, size=4)(2,3)), # x -> - x mirror operation (i.e. yz mirror plane)
                 # dictionary with irrep generators
                               {
                                    # Bradley, Cracknell p. 62
                                    'E': (ra, -lam),
                                    # Bradley, Cracknell p. 59, or Dresselhaus, Table A.18
                                    'A1': (1,1),
                                    'A2': (1,-1),
                                    'B1': (-1,1),
                                    'B2': (-1,-1),    
                               },
                 # function that generates a tuple with svwf representation generators
                               lambda lMax : (qpms.zrotN_tyty(4, lMax), qpms.xflip_tyty(lMax))
    ),
    'D4h' : SVWFPointGroupInfo('D4h',
                 # permutation group generators
                               (Permutation(0,1,2,3, size=6),  # C4 rotation
                                Permutation(0,1, size=6)(2,3), # x -> - x mirror operation (i.e. yz mirror plane)
                                Permutation(4,5, size=6), # horizontal mirror operation z -> -z (i.e. xy mirror plane)
                               ), 
                 # dictionary with irrep generators
                               {    # product of C4v and zflip
                                    "E'": (ra, -lam, epsilon),
                                    "E''":(ra, -lam, -epsilon),
                                    "A1'": (1,1,1),
                                    "A2'": (1,-1,1),
                                    "A1''": (1,-1,-1),
                                    "A2''": (1,1,-1),
                                    "B1'": (-1,1,1),
                                    "B2'": (-1,-1,1),
                                    "B1''": (-1,-1,-1),
                                    "B2''": (-1,1,-1), 
                               },
                 # function that generates a tuple with svwf representation generators
                               lambda lMax : (qpms.zrotN_tyty(4, lMax), qpms.xflip_tyty(lMax), qpms.zflip_tyty(lMax))
    ),
    'D3h' : SVWFPointGroupInfo('D3h',
                 # permutation group generators
                               ( Permutation(0,1,2, size=5),  # C3 rotation
                                 Permutation(0,2, size=5), # vertical mirror
                                 Permutation(3,4, size=5), # horizontal mirror z -> -z (i.e. xy mirror plane)
                               ), 
                 # dictionary with irrep generators
                               {    # Bradley, Cracknell p. 61
                                    "E'" : (alif, lam, epsilon),
                                    "E''" : (alif, lam, -epsilon),
                                    # Bradley, Cracknell p. 59, or Dresselhaus, Table A.14 (p. 482)
                                    "A1'" : (1,1,1),
                                    "A2'" : (1,-1,1),
                                    "A1''" : (1,-1,-1),
                                    "A2''" : (1,1,-1),
                               },
                 # function that generates a tuple with svwf representation generators
                               lambda lMax, vflip: (qpms.zrotN_tyty(3, lMax), qpms.yflip_tyty(lMax) if vflip == 'y' else qpms.xflip_tyty(lMax), qpms.zflip_tyty(lMax))
    ),
}
