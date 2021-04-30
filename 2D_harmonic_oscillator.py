import numpy as np
from qmsolve import Hamiltonian, animate, dynamic_visualize, SingleParticle


#interaction potential
def harmonic_oscillator(particle):

	kx = 100 # measured in eV / (Å**2)
	ky = 100
	return 0.5 * kx * particle.x**2    +    0.5 * ky * particle.y**2



H = Hamiltonian(particles = SingleParticle(), 
				potential = harmonic_oscillator, 
<<<<<<< Updated upstream
				spatial_ndim = 2, N = 200, extent = 15)


energies, eigenstates = H.solve(max_states = 30)

print(energies)
dynamic_visualize(energies, eigenstates)
=======
				spatial_ndim = 2, N = 100, extent = 15)


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
>>>>>>> Stashed changes
