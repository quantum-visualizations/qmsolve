import numpy as np
from qmsolve import Hamiltonian, SingleParticle, init_visualization, Å


#interaction potential
def di_coulomb_potential(particle):
	k_c = 1 # (e*e / (4 * np.pi * epsilon_0)) 

	R = 1.06 * Å

	r1 = np.sqrt((particle.x - R/2)**2 + (particle.y)**2 + (particle.z)**2)
	r1 = np.where(r1 < 0.0001, 0.0001, r1) 

	r2 = np.sqrt((particle.x + R/2)**2 + (particle.y)**2 + (particle.z)**2)
	r2 = np.where(r2 < 0.0001, 0.0001, r2) 


	return - k_c/ r1 - k_c/ r2 + k_c/R 



H = Hamiltonian(particles = SingleParticle(), 
				potential = di_coulomb_potential, 
				spatial_ndim = 3, N = 150, extent = 9*Å)


eigenstates = H.solve(max_states = 5, method ='lobpcg')
print(eigenstates.energies)

visualization = init_visualization(eigenstates)
visualization.plot_eigenstate(0, contrast_vals = [0.01, 0.761])
visualization.animate(contrast_vals = [0.01, 0.761])