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
