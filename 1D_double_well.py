import numpy as np
from qmsolve import Hamiltonian, SingleParticle
from qmsolve import visualize_superpositions, dynamic_visualize


def double_well(particle):
	return 30.0*(particle.x**2 - 1.5**2)**2



H = Hamiltonian(particles = SingleParticle(), 
				potential = double_well, 
				spatial_ndim = 1, N = 300, extent = 20)


energies, eigenstates = H.solve(max_states = 10)

print(energies)
dynamic_visualize(energies, eigenstates)
visualize_superpositions(energies, eigenstates, 10, dt=0.01)