import numpy as np
from qmsolve import Hamiltonian, SingleParticle, init_visualization, Å


#interaction potential
def coulomb_potential(particle):

	k_c = 1.0 # (e*e / (4 * np.pi * epsilon_0))  
	r_spher = np.sqrt((particle.z)**2 + (particle.r)**2)
	return - k_c/ r_spher



H = Hamiltonian(particles = SingleParticle(symmetry = "azimuthal", azimuthal_number_m = 0), 
				potential = coulomb_potential, 
				spatial_ndim = 3, N = 500, extent = 50*Å)


eigenstates = H.solve(max_states = 8, method= 'eigs',)
print(eigenstates.energies)
