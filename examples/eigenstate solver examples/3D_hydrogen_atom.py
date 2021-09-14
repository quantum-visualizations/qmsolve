import numpy as np
from qmsolve import Hamiltonian, SingleParticle, init_visualization, Å, V,m


#interaction potential
def coulomb_potential(particle):

	k_c = 1.0 # (e*e / (4 * np.pi * epsilon_0))  

	r = np.sqrt((particle.x)**2 + (particle.y)**2 + (particle.z)**2)
	r = np.where(r < 0.000001, 0.000001, r)
	external_electric_field = 1e3*V/m #shows Stark effect

	return - k_c/ r + particle.z*external_electric_field





H = Hamiltonian(particles = SingleParticle(), 
				potential = coulomb_potential, 
				spatial_ndim = 3, N = 150, extent = 30*Å)



eigenstates = H.solve( max_states = 14, method ='lobpcg')
print(eigenstates.energies)



visualization = init_visualization(eigenstates)
visualization.plot_eigenstate(5, contrast_vals = [0.01, 0.2])
visualization.animate(contrast_vals = [0.01, 0.2])
