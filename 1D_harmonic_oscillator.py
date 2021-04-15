import numpy as np
from qmsolve import Hamiltonian, animate, dynamic_visualize, SingleParticle


#interaction potential
def harmonic_oscillator(particle):

	k = 100 # measured in eV / (Å**2)
	return 0.5 * k * particle.x**2



H = Hamiltonian(particles = SingleParticle(), 
				potential = harmonic_oscillator, 
				spatial_ndim = 1, N = 300, extent = 20)


energies, eigenstates = H.solve(max_states = 30)

print(energies)
dynamic_visualize(energies, eigenstates)