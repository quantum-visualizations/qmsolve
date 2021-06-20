import numpy as np
from qmsolve import Hamiltonian, SingleParticle, init_visualization


#interaction potential
def coulomb_potential(particle):

	k_c = 14.39964547842567 # (e*e / (4 * np.pi * epsilon_0))  # measured in eV / Å

	r = np.sqrt((particle.x)**2 + (particle.y)**2 + (particle.z)**2)
	r = np.where(r < 0.0001, 0.0001, r)
	external_electric_field = 0.00005 #shows Stark effect

	return - k_c/ r + particle.z*external_electric_field





H = Hamiltonian(particles = SingleParticle(), 
				potential = coulomb_potential, 
				spatial_ndim = 3, N = 150, extent = 40)



eigenstates = H.solve( max_states = 10, N0 = 30, method ='lobpcg-cupy')
print(eigenstates.energies)



visualization = init_visualization(eigenstates)
visualization.plot_eigenstate(5, contrast_vals = [0.01, 0.391])
visualization.animate(contrast_vals = [0.01, 0.391])
