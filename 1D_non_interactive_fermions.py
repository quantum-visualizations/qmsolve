import numpy as np
from qmsolve import Hamiltonian, TwoFermions, init_visualization


def coulomb_interaction(fermions):

	k = 0. # measured in eV * Å**2
	r = (fermions.x1 - fermions.x2)
	r = np.where(r < 0.0001, 0.0001, r)
	V = k/ r**2
	return V

H = Hamiltonian(particles = TwoFermions(), 
				potential = coulomb_interaction, # change this to harmonic_interaction to check what happens!
				spatial_ndim = 1, N = 100, extent = 10)


eigenstates = H.solve(max_states = 90)
print(eigenstates.energies)

visualization = init_visualization(eigenstates)

#visualization.plot_eigenstate(6)
#visualization.slider_plot()
visualization.animate()