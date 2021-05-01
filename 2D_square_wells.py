import numpy as np
from qmsolve import Hamiltonian, SingleParticle, init_visualization


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


eigenstates = H.solve(max_states = 10)
coeffs = np.zeros([10], np.complex128)
coeffs[6] = 1.0
coeffs[7] = 1.0j
visualization = init_visualization(eigenstates)
visualization.superpositions(coeffs, dt=0.02,
							 )
