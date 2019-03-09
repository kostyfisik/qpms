from qpms import Particle, CTMatrix, BaseSpec, FinitePointGroup, ScatteringSystem
from qpms.symmetries import point_group_info
import numpy as np

sym = FinitePointGroup(point_group_info['D3h'])
bspec2 = BaseSpec(lMax=2)
bspec1 = BaseSpec(lMax=1)
t1 = CTMatrix(bspec1, np.diag(np.random.random(len(bspec1))))
t2 = CTMatrix(bspec2, np.diag(np.random.random(len(bspec2))))
p1 = Particle((1,2,),t1)
p2 = Particle((1,2,3),t1)
p3 = Particle((0.1,2),t2)
ss = ScatteringSystem([p1, p2, p3], sym)

#print(ss.fecv_size, ss.saecv_sizes, ss.nirreps, ss.nirrep_names)

fullvector = np.random.rand(ss.fecv_size) + 1j*np.random.rand(ss.fecv_size)
packedvectors = [(iri, ss.pack_vector(fullvector, iri)) for iri in range(ss.nirreps)]
unpackedvectors = np.array([ss.unpack_vector(v[1], v[0]) for v in packedvectors])
rec_fullvector = np.sum(unpackedvectors, axis=0)
thediff = np.amax(abs(rec_fullvector-fullvector))
assert(thediff < 1e-14)

packedmatrices = list()
for iri in range(ss.nirreps):
    d = ss.saecv_sizes[iri]
    m = np.random.rand(d,d)+1j*np.random.rand(d,d)
    packedmatrices.append((iri,m))

fullmatrix = np.zeros((ss.fecv_size, ss.fecv_size), dtype=complex)
for iri, m in packedmatrices:
    fullmatrix += ss.unpack_matrix(m, iri)

for iri, m in packedmatrices:
    print (m.shape)
    repackedmatrix = ss.pack_matrix(fullmatrix,iri)
    print(np.amax(abs(repackedmatrix-m)))
    
