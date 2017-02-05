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
