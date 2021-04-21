import numpy as np
from qmsolve import Hamiltonian,  SingleParticle, init_visualization


#interaction potential
def harmonic_oscillator(particle):

	k = 100 # measured in eV / (Å**2)
	return 0.5 * k * particle.x**2



H = Hamiltonian(particles = SingleParticle(), 
				potential = harmonic_oscillator, 
				spatial_ndim = 1, N = 300, extent = 20)


eigenstates = H.solve(max_states = 30)

print(eigenstates.energies)

visualization = init_visualization(eigenstates)
#visualization.plot_eigenstate(6)
visualization.slider_plot()
#visualization.animate()