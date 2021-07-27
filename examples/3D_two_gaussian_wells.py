import numpy as np
from qmsolve import Hamiltonian, SingleParticle, init_visualization,Å


#interaction potential
def two_gaussian_wells(particle):
	𝜇 = 0.7*Å
	σ = 0.5*Å
	V = 25.72*(2+
	-np.exp((-(particle.x-𝜇)**2 -(particle.y)**2  -(particle.z)**2 ) / (2*σ**2))
	-np.exp((-(particle.x+𝜇)**2 -(particle.y)**2  -(particle.z)**2 ) / (2*σ**2)))
	return V



H = Hamiltonian(particles = SingleParticle(), 
				potential = two_gaussian_wells, 
				spatial_ndim = 3, N = 90, extent = 3*Å)


eigenstates = H.solve( max_states = 50, N0 = 30, method ='lobpcg')
print(eigenstates.energies)

visualization = init_visualization(eigenstates)
visualization.plot_eigenstate(40, contrast_vals = [0.1, 0.25])
visualization.animate(contrast_vals = [0.1, 0.25])
