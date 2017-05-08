import math
import numpy as np
nx = None

_s3 = math.sqrt(3)


def generate_trianglepoints(maxlayer, include_origin = False, v3d = True, circular = True, sixthindices = False, mirrorindices = False):
    _e6 = np.array([[math.cos(2*math.pi*i/6),math.sin(2*math.pi*i/6),0] if v3d else [math.cos(2*math.pi*i/6),math.sin(2*math.pi*i/6)] for i in range(6)])
    points = np.empty((3*maxlayer*(maxlayer+1)+(1 if include_origin else 0), 3 if v3d else 2))
    point_i = 0
    if (include_origin):
        points[0] = np.array((0,0,0) if v3d else (0,0))
        point_i = 1
    if sixthindices:
        si = np.empty((6,(maxlayer*(maxlayer+1))//2), dtype=int)
        sii = [0,0,0,0,0,0]
    if mirrorindices:
        if(maxlayer < 3):
            mi = np.empty((2,0,), dtype=int)
        else:
            #layer indices start from one!
            ilayer = np.arange(1, maxlayer+1) # We need first to "copy" layer indices to correspond to the Muster count
            mustercount = (ilayer - 1)//2
            mustercum = np.cumsum(mustercount)
            layerstart = np.zeros((mustercum[maxlayer - 1]), dtype=int)
            layerstart[mustercum[:(maxlayer-1)]] = 1
            layer = np.cumsum(layerstart) + 2 # That's it
            lb = 3*layer*(layer-1) # layer base (lowest) index
            li = np.arange(len(layer)) - mustercum[layer-2] # muster indices for each layers
            mi = np.empty((2, len(layer)), dtype=int)
            mi[0] = lb + 1 + li + include_origin
            mi[1] = lb + layer - (1 + li) + include_origin
            # there are two non-musters in each even layer, one non-muster in each odd
        layer = (2*np.arange(((3*maxlayer)//2))+1)//3 + 1
        nmi = 3*layer*(layer-1)
        nmi[2::3] += layer[2::3] // 2 # second non-musters in even layers
        nmi += include_origin
        
    for layer in range(1,maxlayer+1):
        for i in range(6):
            base = _e6[i]*layer
            shift = _e6[(i+2)%6]
            ar = np.arange(layer)
            points[point_i:(point_i+layer)] = base[nx,:] + ar[:,nx] * shift[nx,:]
            if sixthindices:
                si[i, sii[i]:sii[i]+layer] = point_i + ar
                sii[i] += layer
            point_i += layer
    if (circular):
        mask = (np.sum(points * points, axis = -1) <= maxlayer * maxlayer * 3/ 4 + 0.1)  # UGLY FIX OF ASYMMETRY BECAUSE OF ROUNDING ERROR      
        points = points[mask]
        if sixthindices:
            cum = np.cumsum(mask) - 1
            mask0 = mask[si[0]]
            si_ = si[:,mask0]
            si = cum[si_]   
        if mirrorindices:
            cum = np.cumsum(mask) - 1
            mask0 = mask[mi[0]]
            mi_ = mi[:,mask0]
            mi = cum[mi_]
            mask0 = mask[nmi]
            nmi_ = nmi[mask0]
            nmi = cum[nmi_]
    if not (mirrorindices or sixthindices):
        return points
    else:
        return {'points': points,
                'si' : si if sixthindices else None,
                'mi' : mi if mirrorindices else None,
                'nmi' : nmi if mirrorindices else None}

def generate_trianglepoints_hexcomplement(maxlayer, v3d = True, circular = True, thirdindices = False, mirrorindices=False):
    _e6 = np.array([[math.cos(2*math.pi*i/6),math.sin(2*math.pi*i/6),0] if v3d else [math.cos(2*math.pi*i/6),math.sin(2*math.pi*i/6)] for i in range(6)])
    _f6 = np.array([[-math.sin(2*math.pi*i/6),math.cos(2*math.pi*i/6),0] if v3d else [math.sin(2*math.pi*i/6),-math.cos(2*math.pi*i/6)] for i in range(6)])
    points = np.empty((3*maxlayer*maxlayer, 3 if v3d else 2))
    point_i = 0
    # 3 * layer ** 2 is the basis index for a layer, a layer contains 3 * (2*layer + 1) points
    if thirdindices:
        ii = np.arange(maxlayer**2)
        layer = np.empty((maxlayer**2), dtype=int)
        layer = np.sqrt(ii, out=layer, casting='unsafe')
        #ti0 = 2*layer**2 + ii 
        ti = np.arange(3)[:, nx] * (2*layer + 1)[nx, :] + (2*layer**2 + ii)[nx,:]
    if mirrorindices:
        ii = np.arange(((maxlayer-1)*maxlayer)//2)
        layer = np.empty((((maxlayer-1)*maxlayer)//2), dtype=int)
        layer = (np.sqrt(1+8*ii, out=layer, casting='unsafe')+1)//2
        li = ii - ((layer ) * (layer-1))//2# numbers indices in a each layer
        lb = 3*layer **2  # base index of a layer
        mi = np.empty((2,len(layer)), dtype=int)
        mi[0] = lb + li + layer % 2
        mi[1] = lb + 2*layer - li
        # indices of non-mirrored/self-mirrored
        layer = np.arange(maxlayer)
        lb = 3 * layer**2
        nmi = lb + ((layer + 1) % 2) * layer
    for layer in range(0,maxlayer):
        if (layer % 2): # odd layer
            for i in range(3):
                base = _f6[(2*i-1)%6] * ((0.5 + 1.5 * layer) / _s3)
                shift = _e6[(2*i+2)%6]
                count = (layer + 1) // 2
                ar = np.arange(count)
                points[point_i:point_i+count,:] = base + ar[:,nx]*shift[nx,:]
                point_i += count
                
                base = _e6[(2*i+1)%6]*layer + _f6[(2*i)%6] / _s3
                shift = _e6[(2*i+3)%6]
                count = layer
                ar = np.arange(count)
                points[point_i:point_i+count,:] = base + ar[:,nx]*shift[nx,:]
                point_i += count

                base = _e6[(2*i+2)%6]*layer + _f6[(2*i)%6] / _s3
                shift = _e6[(2*i+4)%6]
                count = (layer + 1) // 2
                ar = np.arange(count)
                points[point_i:point_i+count,:] = base + ar[:,nx]*shift[nx,:]
                point_i += count
        else: # even layer:
            for i in range(3):
                shift = _e6[(2*i+2)%6]
                base = _f6[(2*i-1)%6] * ((0.5 + 1.5 * layer) / _s3) + shift / 2
                count = layer // 2
                ar = np.arange(count)
                points[point_i:point_i+count,:] = base + ar[:,nx]*shift[nx,:]
                point_i += count
                
                base = _e6[(2*i+1)%6]*layer + _f6[(2*i)%6] / _s3
                shift = _e6[(2*i+3)%6]
                count = layer
                ar = np.arange(count)
                points[point_i:point_i+count,:] = base + ar[:,nx]*shift[nx,:]
                point_i += count
                
                base = _e6[(2*i+2)%6]*layer + _f6[(2*i)%6] / _s3
                shift = _e6[(2*i+4)%6]
                count = (layer + 2) // 2
                ar = np.arange(count)
                points[point_i:point_i+count,:] = base + ar[:,nx]*shift[nx,:]
                point_i += count
    #if (mirrorindices):
        
    if (circular):
        mask = (np.sum(points * points, axis = -1) <= maxlayer * maxlayer * 3/ 4 + 0.01)  # UGLY FIX OF ASYMMETRY BECAUSE OF ROUNDING ERROR      
        points = points[mask]
        if thirdindices:
            cum = np.cumsum(mask) - 1
            mask0 = mask[ti[0]]
            ti_ = ti[:,mask0]
            ti = cum[ti_]  
        if mirrorindices:
            cum = np.cumsum(mask) - 1
            mask0 = mask[mi[0]]
            mi_ = mi[:,mask0]
            mi = cum[mi_]
            mask0 = mask[nmi]
            nmi_ = nmi[mask0]
            nmi = cum[nmi_]
    if not (mirrorindices or thirdindices):
        return points
    else:
        return {'points': points,
                'ti' : ti if thirdindices else None,
                'mi' : mi if mirrorindices else None,
                'nmi' : nmi if mirrorindices else None
               }


from qpms_c import get_mn_y, trans_calculator
from .qpms_p import cart2sph

def hexlattice_precalc_AB_save(file, lMax, k_hexside, maxlayer, circular=True, savepointinfo = False, J_scat=3):
    params = {
        'lMax' : lMax,
        'k_hexside' : k_hexside,
        'maxlayer' : maxlayer,
        'circular' : circular,
        'savepointinfo' : savepointinfo,
        'J_scat' : J_scat
    }
    tpdict = generate_trianglepoints(maxlayer, v3d=True, circular=circular, sixthindices=True, mirrorindices=True)
    tphcdict = generate_trianglepoints_hexcomplement(maxlayer, v3d=True, circular=circular, thirdindices=True, mirrorindices=True)
    my, ny = get_mn_y(lMax)
    nelem = len(my)
    a_self_nm = np.empty((tpdict['nmi'].shape[0],nelem,nelem), dtype=complex)
    b_self_nm = np.empty((tpdict['nmi'].shape[0],nelem,nelem), dtype=complex)
    a_self_m0 = np.empty((tpdict['mi'].shape[1],nelem,nelem), dtype=complex)
    b_self_m0 = np.empty((tpdict['mi'].shape[1],nelem,nelem), dtype=complex)
    a_d2u_nm = np.empty((tphcdict['nmi'].shape[0],nelem,nelem), dtype=complex)
    b_d2u_nm = np.empty((tphcdict['nmi'].shape[0],nelem,nelem), dtype=complex)
    a_d2u_m0 = np.empty((tphcdict['mi'].shape[1],nelem,nelem), dtype=complex)
    b_d2u_m0 = np.empty((tphcdict['mi'].shape[1],nelem,nelem), dtype=complex)
    
    k_0 = k_hexside*_s3 # not really a wave vector here because of the normalisation!
    tc = trans_calculator(lMax)
    
    y = np.arange(nelem)

    points = tpdict['points'][tpdict['nmi']]
    d_i2j = cart2sph(points)
    a_self_nm, b_self_nm = tc.get_AB_arrays(k_0*d_i2j[:,0],d_i2j[:,1],d_i2j[:,2],np.array([False]),J_scat)

    points = tpdict['points'][tpdict['mi'][0]]
    d_i2j = cart2sph(points)
    a_self_m0, b_self_m0 = tc.get_AB_arrays(k_0*d_i2j[:,0],d_i2j[:,1],d_i2j[:,2],np.array([False]),J_scat)
   
    points = tphcdict['points'][tphcdict['nmi']]
    d_i2j = cart2sph(points)
    a_d2u_nm, b_d2u_nm = tc.get_AB_arrays(k_0*d_i2j[:,0],d_i2j[:,1],d_i2j[:,2],np.array([False]),J_scat)
  
    points = tphcdict['points'][tphcdict['mi'][0]]
    d_i2j = cart2sph(points)
    a_d2u_m0, b_d2u_m0 = tc.get_AB_arrays(k_0*d_i2j[:,0],d_i2j[:,1],d_i2j[:,2],np.array([False]),J_scat)
  
    tosave = {
        'a_self_nm' : a_self_nm,
        'a_self_m0' : a_self_m0,
        'b_self_nm' : b_self_nm,
        'b_self_m0' : b_self_m0,
        'a_d2u_nm' : a_d2u_nm,
        'a_d2u_m0' : a_d2u_m0,
        'b_d2u_nm' : b_d2u_nm,
        'b_d2u_m0' : b_d2u_m0,
        'precalc_params' : params
    }
    if savepointinfo:
        tosave['tp_points'] = tpdict['points'],
        tosave['tp_si'] = tpdict['si'],
        tosave['tp_mi'] = tpdict['mi'],
        tosave['tp_nmi'] = tpdict['nmi']
        tosave['tphc_points'] = tphcdict['points'],
        tosave['tphc_ti'] = tphcdict['ti'],
        tosave['tphc_mi'] = tphcdict['mi'],
        tosave['tphc_nmi'] = tphcdict['nmi']
    np.savez(file, **tosave)

def hexlattice_precalc_AB_loadunwrap(file, tpdict = None, tphcdict = None, return_points = False):
    npz = np.load(file)
    precalc_params = npz['precalc_params'][()]
    my, ny = get_mn_y(precalc_params['lMax'])
    nelem = len(my)
    # this I should have made more universal...
    if precalc_params['savepointinfo']:
        if not tpdict:
            tpdict = {
                'points' : npz['tp_points'],
                'si' : npz['tp_si'],
                'mi' : npz['tp_mi'],
                'nmi' : npz['tp_nmi'],
            }
        if not tphcdict:
            tphcdict = {
                'points' : npz['tphc_points'],
                'ti' : npz['tphc_ti'],
                'mi' : npz['tphc_mi'],
                'nmi' : npz['tphc_nmi']            
        }
    else:
        if not tpdict:
            tpdict = generate_trianglepoints(maxlayer = precalc_params['maxlayer'], v3d=True, 
                                             circular=precalc_params['circular'], sixthindices=True, mirrorindices=True)
        if not tphcdict:
            tphcdict = generate_trianglepoints_hexcomplement(maxlayer=precalc_params['maxlayer'], v3d=True, 
                                                             circular=precalc_params['circular'], thirdindices=True, mirrorindices=True)
        
    # For some obscure reason, I keep getting trailing single-dimension in the beginning for these arrays
    for a in (tpdict['points'], tphcdict['points'], tpdict['si'], tpdict['mi'],
             tphcdict['ti'], tphcdict['mi']):
        if len(a.shape) > 2:
            a.shape = a.shape[1::] 
        
    self_tr = tpdict['points'] 
    d2u_tr = tphcdict['points']
    if len(self_tr.shape)>2:
        self_tr = np.reshape(self_tr, self_tr.shape[1::])
    if len(d2u_tr.shape)>2:
        d2u_tr = np.reshape(d2u_tr, d2u_tr.shape[1::])
    u2d_tr = -d2u_tr
    a_self = np.empty((self_tr.shape[0],nelem,nelem), dtype=complex)
    b_self = np.empty((self_tr.shape[0],nelem,nelem), dtype=complex)
    a_d2u  = np.empty(( d2u_tr.shape[0],nelem,nelem), dtype=complex)
    b_d2u  = np.empty(( d2u_tr.shape[0],nelem,nelem), dtype=complex)
    a_self[tpdict['nmi']]=npz['a_self_nm']
    a_self[tpdict['mi'][0]]=npz['a_self_m0']
    b_self[tpdict['nmi']]=npz['b_self_nm']
    b_self[tpdict['mi'][0]]=npz['b_self_m0']
    mirrorangles = cart2sph(self_tr[tpdict['mi'][1]])[:,2] - cart2sph(self_tr[tpdict['mi'][0]])[:,2]
    a_self[tpdict['mi'][1],:,:] = a_self[tpdict['mi'][0],:,:] * np.exp(1j*mirrorangles[:,nx,nx]*(my[nx,nx,:]-my[nx,:,nx]))
    b_self[tpdict['mi'][1],:,:] = b_self[tpdict['mi'][0],:,:] * np.exp(1j*mirrorangles[:,nx,nx]*(my[nx,nx,:]-my[nx,:,nx]))
    for i in range(1,6):
        a_self[tpdict['si'][i],:,:] = a_self[tpdict['si'][0],:,:] * np.exp(1j*math.pi/3*i*(my[nx,:]-my[:,nx]))
        b_self[tpdict['si'][i],:,:] = b_self[tpdict['si'][0],:,:] * np.exp(1j*math.pi/3*i*(my[nx,:]-my[:,nx]))
    a_d2u[tphcdict['nmi']]=npz['a_d2u_nm']
    a_d2u[tphcdict['mi'][0]]=npz['a_d2u_m0']
    b_d2u[tphcdict['nmi']]=npz['b_d2u_nm']
    b_d2u[tphcdict['mi'][0]]=npz['b_d2u_m0']
    mirrorangles = cart2sph(self_tr[tphcdict['mi'][1]])[:,2] - cart2sph(self_tr[tphcdict['mi'][0]])[:,2]
    a_d2u[tphcdict['mi'][1],:,:] = a_d2u[tphcdict['mi'][0],:,:] * np.exp(1j*mirrorangles[:,nx,nx]*(my[nx,nx,:]-my[nx,:,nx]))
    b_d2u[tphcdict['mi'][1],:,:] = b_d2u[tphcdict['mi'][0],:,:] * np.exp(1j*mirrorangles[:,nx,nx]*(my[nx,nx,:]-my[nx,:,nx]))
    for i in (1,-1):
        a_d2u[tphcdict['ti'][i],:,:] = a_d2u[tphcdict['ti'][0],:,:] * np.exp(i*2j*math.pi/3*(my[nx,:]-my[:,nx]))
        b_d2u[tphcdict['ti'][i],:,:] = b_d2u[tphcdict['ti'][0],:,:] * np.exp(i*2j*math.pi/3*(my[nx,:]-my[:,nx]))
    a_u2d = a_d2u * (-1)**(my[nx,:]-my[:,nx])
    b_u2d = b_d2u * (-1)**(my[nx,:]-my[:,nx])
    d = {
        'a_self' : a_self,
        'b_self' : b_self,
        'a_d2u'  : a_d2u,
        'b_d2u'  : b_d2u,
        'a_u2d'  : a_u2d,
        'b_u2d'  : b_u2d,
    }
    for k in precalc_params.keys():
        d[k] = precalc_params[k]
    if return_points:
        d['d2u_tr'] = tphcdict['points']
        d['u2d_tr'] = -tphcdict['points']
        d['self_tr'] = tpdict['points']
    return d
    
def hexlattice_get_AB(lMax, k_hexside, maxlayer, circular=True, return_points = True, J_scat=3):
    params = {
        'lMax' : lMax,
        'k_hexside' : k_hexside,
        'maxlayer' : maxlayer,
        'circular' : circular,
        'savepointinfo' : return_points, # should I delete this key?
        'J_scat' : J_scat
    }
    tpdict = generate_trianglepoints(maxlayer, v3d=True, circular=circular, sixthindices=True, mirrorindices=True)
    tphcdict = generate_trianglepoints_hexcomplement(maxlayer, v3d=True, circular=circular, thirdindices=True, mirrorindices=True)
    my, ny = get_mn_y(lMax)
    nelem = len(my)
    a_self_nm = np.empty((tpdict['nmi'].shape[0],nelem,nelem), dtype=complex)
    b_self_nm = np.empty((tpdict['nmi'].shape[0],nelem,nelem), dtype=complex)
    a_self_m0 = np.empty((tpdict['mi'].shape[1],nelem,nelem), dtype=complex)
    b_self_m0 = np.empty((tpdict['mi'].shape[1],nelem,nelem), dtype=complex)
    a_d2u_nm = np.empty((tphcdict['nmi'].shape[0],nelem,nelem), dtype=complex)
    b_d2u_nm = np.empty((tphcdict['nmi'].shape[0],nelem,nelem), dtype=complex)
    a_d2u_m0 = np.empty((tphcdict['mi'].shape[1],nelem,nelem), dtype=complex)
    b_d2u_m0 = np.empty((tphcdict['mi'].shape[1],nelem,nelem), dtype=complex)
    
    k_0 = k_hexside*_s3 # not really a wave vector here because of the normalisation!
    tc = trans_calculator(lMax)
    
    y = np.arange(nelem)

    points = tpdict['points'][tpdict['nmi']]
    d_i2j = cart2sph(points)
    a_self_nm, b_self_nm = tc.get_AB_arrays(k_0*d_i2j[:,0],d_i2j[:,1],d_i2j[:,2],np.array([False]),J_scat)

    points = tpdict['points'][tpdict['mi'][0]]
    d_i2j = cart2sph(points)
    a_self_m0, b_self_m0 = tc.get_AB_arrays(k_0*d_i2j[:,0],d_i2j[:,1],d_i2j[:,2],np.array([False]),J_scat)
   
    points = tphcdict['points'][tphcdict['nmi']]
    d_i2j = cart2sph(points)
    a_d2u_nm, b_d2u_nm = tc.get_AB_arrays(k_0*d_i2j[:,0],d_i2j[:,1],d_i2j[:,2],np.array([False]),J_scat)
  
    points = tphcdict['points'][tphcdict['mi'][0]]
    d_i2j = cart2sph(points)
    a_d2u_m0, b_d2u_m0 = tc.get_AB_arrays(k_0*d_i2j[:,0],d_i2j[:,1],d_i2j[:,2],np.array([False]),J_scat)
    '''
    tosave = {
        'a_self_nm' : a_self_nm,
        'a_self_m0' : a_self_m0,
        'b_self_nm' : b_self_nm,
        'b_self_m0' : b_self_m0,
        'a_d2u_nm' : a_d2u_nm,
        'a_d2u_m0' : a_d2u_m0,
        'b_d2u_nm' : b_d2u_nm,
        'b_d2u_m0' : b_d2u_m0,
        'precalc_params' : params
    }
    if savepointinfo:
        tosave['tp_points'] = tpdict['points'],
        tosave['tp_si'] = tpdict['si'],
        tosave['tp_mi'] = tpdict['mi'],
        tosave['tp_nmi'] = tpdict['nmi']
        tosave['tphc_points'] = tphcdict['points'],
        tosave['tphc_ti'] = tphcdict['ti'],
        tosave['tphc_mi'] = tphcdict['mi'],
        tosave['tphc_nmi'] = tphcdict['nmi']
    np.savez(file, **tosave)
    '''
    self_tr = tpdict['points'] 
    d2u_tr = tphcdict['points']
    if len(self_tr.shape)>2:
        self_tr = np.reshape(self_tr, self_tr.shape[1::])
    if len(d2u_tr.shape)>2:
        d2u_tr = np.reshape(d2u_tr, d2u_tr.shape[1::])
    u2d_tr = -d2u_tr
    a_self = np.empty((self_tr.shape[0],nelem,nelem), dtype=complex)
    b_self = np.empty((self_tr.shape[0],nelem,nelem), dtype=complex)
    a_d2u  = np.empty(( d2u_tr.shape[0],nelem,nelem), dtype=complex)
    b_d2u  = np.empty(( d2u_tr.shape[0],nelem,nelem), dtype=complex)
    a_self[tpdict['nmi']]=a_self_nm
    a_self[tpdict['mi'][0]]=a_self_m0
    b_self[tpdict['nmi']]=b_self_nm
    b_self[tpdict['mi'][0]]=b_self_m0
    mirrorangles = cart2sph(self_tr[tpdict['mi'][1]])[:,2] - cart2sph(self_tr[tpdict['mi'][0]])[:,2]
    a_self[tpdict['mi'][1],:,:] = a_self[tpdict['mi'][0],:,:] * np.exp(1j*mirrorangles[:,nx,nx]*(my[nx,nx,:]-my[nx,:,nx]))
    b_self[tpdict['mi'][1],:,:] = b_self[tpdict['mi'][0],:,:] * np.exp(1j*mirrorangles[:,nx,nx]*(my[nx,nx,:]-my[nx,:,nx]))
    for i in range(1,6):
        a_self[tpdict['si'][i],:,:] = a_self[tpdict['si'][0],:,:] * np.exp(1j*math.pi/3*i*(my[nx,:]-my[:,nx]))
        b_self[tpdict['si'][i],:,:] = b_self[tpdict['si'][0],:,:] * np.exp(1j*math.pi/3*i*(my[nx,:]-my[:,nx]))
    a_d2u[tphcdict['nmi']]=a_d2u_nm
    a_d2u[tphcdict['mi'][0]]=a_d2u_m0
    b_d2u[tphcdict['nmi']]=b_d2u_nm
    b_d2u[tphcdict['mi'][0]]=b_d2u_m0
    mirrorangles = cart2sph(self_tr[tphcdict['mi'][1]])[:,2] - cart2sph(self_tr[tphcdict['mi'][0]])[:,2]
    a_d2u[tphcdict['mi'][1],:,:] = a_d2u[tphcdict['mi'][0],:,:] * np.exp(1j*mirrorangles[:,nx,nx]*(my[nx,nx,:]-my[nx,:,nx]))
    b_d2u[tphcdict['mi'][1],:,:] = b_d2u[tphcdict['mi'][0],:,:] * np.exp(1j*mirrorangles[:,nx,nx]*(my[nx,nx,:]-my[nx,:,nx]))
    for i in (1,-1):
        a_d2u[tphcdict['ti'][i],:,:] = a_d2u[tphcdict['ti'][0],:,:] * np.exp(i*2j*math.pi/3*(my[nx,:]-my[:,nx]))
        b_d2u[tphcdict['ti'][i],:,:] = b_d2u[tphcdict['ti'][0],:,:] * np.exp(i*2j*math.pi/3*(my[nx,:]-my[:,nx]))
    a_u2d = a_d2u * (-1)**(my[nx,:]-my[:,nx])
    b_u2d = b_d2u * (-1)**(my[nx,:]-my[:,nx])
    d = {
        'a_self' : a_self,
        'b_self' : b_self,
        'a_d2u'  : a_d2u,
        'b_d2u'  : b_d2u,
        'a_u2d'  : a_u2d,
        'b_u2d'  : b_u2d,
    }
    for k in params.keys():
        d[k] = params[k]
    if return_points:
        d['d2u_tr'] = tphcdict['points']
        d['u2d_tr'] = -tphcdict['points']
        d['self_tr'] = tpdict['points']
    return d
 

