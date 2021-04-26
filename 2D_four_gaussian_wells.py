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
visualization.animate()

