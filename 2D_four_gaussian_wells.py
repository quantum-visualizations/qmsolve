import numpy as np
from qmsolve import Hamiltonian, SingleParticle, init_visualization


#interaction potential
def four_gaussian_wells(particle):
	𝜇 = 1.6
	σ = 0.5
	V = 600*(4-np.exp((-(particle.x)**2 -(particle.y-𝜇)**2 ) / (2*σ**2))
	-np.exp((-(particle.x-𝜇)**2 -(particle.y)**2 ) / (2*σ**2))
	-np.exp((-(particle.x+𝜇)**2 -(particle.y)**2 ) / (2*σ**2))
	-np.exp((-(particle.x)**2 -(particle.y+𝜇)**2 ) / (2*σ**2)))
	return V



H = Hamiltonian(particles = SingleParticle(), 
				potential = four_gaussian_wells, 
				spatial_ndim = 2, N = 100, extent = 8)


eigenstates = H.solve(max_states = 60)

print(eigenstates.energies)
visualization = init_visualization(eigenstates)
#visualization.plot_eigenstate(6)
coeffs = np.zeros([30], np.complex128)
coeffs[6] = 1.0
coeffs[7] = 1.0j
coeffs[28] = (1.0 - 1.0j)/np.sqrt(2.0)
coeffs[29] = (1.0 + 1.0j)/np.sqrt(2.0)
visualization.superpositions(coeffs, dt=0.0002, 
							 # hide_controls=True, 
							 # save_animation=True, frames=7*30
							 )
# visualization.animate()

