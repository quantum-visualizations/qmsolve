import numpy as np
from qmsolve import Hamiltonian, SingleParticle, init_visualization



#interaction potential
def harmonic_oscillator(particle):

	kx = 2 # measured in eV / (Å**2)
	ky = 2
	return 0.5 * kx * particle.x**2    +    0.5 * ky * particle.y**2



H = Hamiltonian(particles = SingleParticle(), 
				potential = harmonic_oscillator, 
				spatial_ndim = 2, N = 200, extent = 15)


eigenstates = H.solve(max_states = 30)

print(eigenstates.energies)
visualization = init_visualization(eigenstates)
#visualization.plot_eigenstate(6)
#visualization.slider_plot()
#visualization.animate()
visualization.superpositions(10, dt=0.04)