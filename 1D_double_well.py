import numpy as np
from qmsolve import Hamiltonian,  SingleParticle, init_visualization



def double_well(particle):
	return 30.0*(particle.x**2 - 1.5**2)**2



H = Hamiltonian(particles = SingleParticle(), 
				potential = double_well, 
				spatial_ndim = 1, N = 300, extent = 20)


eigenstates = H.solve(max_states = 10)

vis = init_visualization(eigenstates)
coeffs = np.zeros([10], np.complex128)
coeffs[0], coeffs[1] = 1.0, 1.0j
vis.superpositions(states=coeffs, dt=1.0)