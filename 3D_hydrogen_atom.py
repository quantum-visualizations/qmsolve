import numpy as np
from qmsolve import Hamiltonian, SingleParticle, init_visualization


#interaction potential
def coulomb_potential(particle):

	k_c = 14.39964547842567 # (e*e / (4 * np.pi * epsilon_0))  # measured in eV / Å

	r = np.sqrt((particle.x)**2 + (particle.y)**2 + (particle.z)**2)
	r = np.where(r < 0.0001, 0.0001, r)
	return - k_c/ r




H = Hamiltonian(particles = SingleParticle(), 
				potential = coulomb_potential, 
				spatial_ndim = 3, N = 30, extent = 30)


eigenstates = H.solve(max_states = 20)
print(eigenstates.energies)

visualization = init_visualization(eigenstates)
visualization.plot_eigenstate(19, contrast_vals = [0.04, 0.09])
visualization.animate(contrast_vals = [0.04, 0.09])