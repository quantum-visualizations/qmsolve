import numpy as np
from qmsolve import Hamiltonian, SingleParticle, init_visualization,Å



#interaction potential
def harmonic_oscillator(particle):

	kx = 0.02 
	ky = 0.02
	return 0.5 * kx * particle.x**2    +    0.5 * ky * particle.y**2



H = Hamiltonian(particles = SingleParticle(), 
				potential = harmonic_oscillator, 
				spatial_ndim = 2, N = 100, extent = 15*Å)


eigenstates = H.solve(max_states = 30)

print(eigenstates.energies)
visualization = init_visualization(eigenstates)
#visualization.plot_eigenstate(6)
visualization.slider_plot()
#visualization.animate()
coeffs = np.zeros([10], np.complex128)
coeffs[1] = 1.0
coeffs[2] = 1.0j
visualization.superpositions(coeffs, dt=0.01, 
							 # save_animation=True, frames=60
							 )