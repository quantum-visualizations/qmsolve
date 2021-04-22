import numpy as np
from qmsolve import Hamiltonian, TwoFermions, init_visualization


def harmonic_interaction(fermions):

	k = 50 # measured in eV / (Å**2)

	l0 = 5 # measured in Å
	V = 0.5*k*(fermions.x1 - fermions.x2 - l0) **2
	return V



def coulomb_interaction(fermions):

	k = 500. # measured in eV * Å**2
	r = (fermions.x1 - fermions.x2)
	r = np.where(r < 0.0001, 0.0001, r)
	V = k/ r**2
	return V

H = Hamiltonian(particles = TwoFermions(), 
				potential = coulomb_interaction, # change this to harmonic_interaction to check what happens!
				spatial_ndim = 1, N = 200, extent = 10)


eigenstates = H.solve(max_states = 90)
print(eigenstates.energies)

visualization = init_visualization(eigenstates)

#visualization.slider_plot()
visualization.animate()