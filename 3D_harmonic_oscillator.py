import numpy as np
from qmsolve import visualization
from qmsolve import Hamiltonian, SingleParticle, Eigenstates


#interaction potential
def harmonic_oscillator(particle):

	kx, ky, kz = [2]*3 # measured in eV / (Å**2)
	return 0.500 * kx * particle.x**2 +\
		   0.51 * ky * particle.y**2 +\
		   0.52 * kz * particle.z**2


H = Hamiltonian(particles=SingleParticle(),
                  N=30, potential=harmonic_oscillator, 
                  extent=10, spatial_ndim=3)


eigenstates = H.solve(max_states = 32)
print(eigenstates.energies)

visualization = visualization.init_visualization(eigenstates)
visualization.plot_eigenstate(26, contrast_vals = [0.001, 0.5])
#visualization.animate()
