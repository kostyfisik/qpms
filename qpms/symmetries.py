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
    
epsilon = np.eye(2)
alif = np.array(((-1/2,-sqrt(3)/2),(sqrt(3)/2,-1/2)))
bih = np.array(((-1/2,sqrt(3)/2),(-sqrt(3)/2,-1/2)))
lam = np.array(((1,0),(0,-1)))
mim =  np.array(((-1/2,-sqrt(3)/2),(-sqrt(3)/2,1/2)))
nun =  np.array(((-1/2,sqrt(3)/2),(sqrt(3)/2,1/2)))


# Group D3h
# Note that the size argument of permutations is necessary, otherwise e.g. c*c and  b*b would not be evaluated equal
# N.B. the weird elements as Permutation(N) – it means identity permutation of size N+1.
rot3_perm = Permutation(0,1,2, size=5) # C3 rotation
xflip_perm = Permutation(0,2, size=5) # vertical mirror
zflip_perm = Permutation(3,4, size=5) # horizontal mirror
D3h_srcgens = [rot3_perm,xflip_perm,zflip_perm]
D3h_permgroup = PermutationGroup(rot3_perm,xflip_perm,zflip_perm) # D3h

#srcgens = [a,b,c]
D3h_irreps = {
    # Bradley, Cracknell p. 61
    'E1' : generate_grouprep(D3h_permgroup, epsilon, D3h_srcgens, [alif, lam, epsilon], immultop = np.dot, imcmp = np.allclose),
    'E2' : generate_grouprep(D3h_permgroup, epsilon, D3h_srcgens, [alif, lam, -epsilon], immultop = np.dot, imcmp = np.allclose),
    # Bradley, Cracknell p. 59,
    'A1p' : generate_grouprep(D3h_permgroup, 1, D3h_srcgens, [1,1,1]),
    'A2p' : generate_grouprep(D3h_permgroup, 1, D3h_srcgens, [1,-1,1]),
    'A1pp' : generate_grouprep(D3h_permgroup, 1, D3h_srcgens, [1,1,-1]),
    'A2pp' : generate_grouprep(D3h_permgroup, 1, D3h_srcgens, [1,-1,-1]),
}


def mmult_tyty(a, b):
        return(qpms.apply_ndmatrix_left(a, b, (-4,-3)))
def mmult_ptypty(a, b):
    return(qpms.apply_ndmatrix_left(a, b, (-6,-5,-4)))
    
#TODO lepší název fce
def gen_point_D3h_svwf_rep(lMax):
    my, ny = qpms.get_mn_y(lMax)
    nelem = len(my)
    C3_yy = qpms.WignerD_yy_fromvector(lMax, np.array([0,0,2*pi/3]))
    C3_tyty = np.moveaxis(np.eye(2)[:,:,ň,ň] * C3_yy, 2,1)
    zfl_tyty = qpms.zflip_tyty(lMax)
    yfl_tyty = qpms.yflip_tyty(lMax)
    xfl_tyty = qpms.xflip_tyty(lMax)
    I_tyty = np.moveaxis(np.eye(2)[:,:,ň,ň] * np.eye(nelem), 2,1)
    order = D3h_permgroup.order()
    sphrep_full = generate_grouprep(D3h_permgroup, I_tyty, D3h_srcgens, [C3_tyty, xfl_tyty, zfl_tyty], 
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
        
def gen_hexlattice_Kpoint_svwf_rep(lMax, psi):
    my, ny = qpms.get_mn_y(lMax)
    nelem = len(my)
    C3_yy = qpms.WignerD_yy_fromvector(lMax, np.array([0,0,2*pi/3]))
    C3_tyty = np.moveaxis(np.eye(2)[:,:,ň,ň] * C3_yy, 2,1)
    zfl_tyty = qpms.zflip_tyty(lMax)
    yfl_tyty = qpms.yflip_tyty(lMax)
    xfl_tyty = qpms.xflip_tyty(lMax)
    I_tyty = np.moveaxis(np.eye(2)[:,:,ň,ň] * np.eye(nelem), 2,1)
    hex_C3_K_ptypty = np.diag([exp(-psi*1j*2*pi/3),exp(+psi*1j*2*pi/3)])[:,ň,ň,:,ň,ň] * C3_tyty[ň,:,:,ň,:,:]
    hex_zfl_ptypty = np.eye(2)[:,ň,ň,:,ň,ň] * zfl_tyty[ň,:,:,ň,:,:]
    hex_xfl_ptypty = np.array([[0,1],[1,0]])[:,ň,ň,:,ň,ň] * xfl_tyty[ň,:,:,ň,:,:]
    hex_I_ptypty = np.eye((2*2*nelem)).reshape((2,2,nelem,2,2,nelem))
    order = D3h_permgroup.order()
    hex_K_sphrep_full = generate_grouprep(D3h_permgroup, hex_I_ptypty, D3h_srcgens, [hex_C3_K_ptypty, hex_xfl_ptypty, hex_zfl_ptypty], 
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

def gen_hexlattice_Kpoint_svwf_rep_projectors(lMax,psi):
    nelem = lMax * (lMax+2)
    projectors = dict()
    for repi, W in gen_hexlattice_Kpoint_svwf_rep(lMax,psi).items():
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
        for v in tmplist:
            theprojector += (v[:,:,:,ň,ň,ň] * v.conjugate()[ň,ň,ň,:,:,:]).real # TODO check is it possible to have imaginary elements?
        for x in [0, 1, -1,sqrt(0.5),-sqrt(0.5),0.5,-0.5]:
            theprojector[np.isclose(theprojector,x)]=x
        projectors[repi] = theprojector
    return projectors

