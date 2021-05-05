import numpy as np
from qmsolve import Hamiltonian, SingleParticle, init_visualization


#interaction potential
def two_gaussian_wells(particle):
	𝜇 = 0.7
	σ = 0.5
	V = 700*(2+
	-np.exp((-(particle.x-𝜇)**2 -(particle.y)**2  -(particle.z)**2 ) / (2*σ**2))
	-np.exp((-(particle.x+𝜇)**2 -(particle.y)**2  -(particle.z)**2 ) / (2*σ**2)))
	return V



H = Hamiltonian(particles = SingleParticle(), 
				potential = two_gaussian_wells, 
				spatial_ndim = 3, N = 30, extent = 3)


eigenstates = H.solve(max_states = 50)
print(eigenstates.energies)

visualization = init_visualization(eigenstates)
visualization.plot_eigenstate(26)
visualization.animate()