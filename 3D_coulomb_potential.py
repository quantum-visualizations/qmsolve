import numpy as np
from qmsolve import Hamiltonian, SingleParticle
from qmsolve import init_visualization
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
e = H.solve(max_states = 10)
visualization = init_visualization(e)
visualization.animate()
