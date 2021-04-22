import numpy as np
from qmsolve import Hamiltonian, SingleParticle
from qmsolve import animate, dynamic_visualize, visualize
from qmsolve import visualize_superpositions


def double_well(particle):
    N = len(particle.x)
    V = np.zeros([N, N])
    height = 150.0
    V[:, 64*N//128: 65*N//128] = height
    V[64*N//128: 65*N//128, :] = height
    V[0: 4*N//32, :] = height
    V[:, 0: 4*N//32] = height
    V[N - 4*N//32: N, :] = height
    V[:, N - 4*N//32: N] = height
    return V



H = Hamiltonian(particles = SingleParticle(), 
				potential = double_well, 
				spatial_ndim = 2, N = 200, extent = 20)


energies, eigenstates = H.solve(max_states = 10)

print(energies)
visualize_superpositions(energies, eigenstates, 10, dt=0.1)