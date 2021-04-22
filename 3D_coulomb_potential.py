import numpy as np
from qmsolve import Hamiltonian, SingleParticle
from qmsolve import visualize3D, animate3D
from qmsolve import k_c


spatial_ndim = 1
L = 10.0
N = 32


def coulomb(particle):
    r = particle.x, particle.y, particle.z
    return -k_c/np.sqrt(sum([r_i**2 for r_i in r]))


H = Hamiltonian(particles=SingleParticle(),
                N=N, potential=coulomb, 
                extent=L, spatial_ndim=3)
e, v = H.solve(max_states = 10)


print(e)
visualize3D(e, v, 9)
visualize3D(e, v, 9, plot_type='contour')
animate3D(e, v)
