import numpy as np
from qmsolve import Hamiltonian, TwoFermions, init_visualization,Å


def harmonic_interaction(fermions):

	k = 0.5

	l0 = 5 # measured in Å
	V = 0.5*k*(fermions.x1 - fermions.x2 - l0) **2
	return V



def coulomb_interaction(fermions):

	k = 34.
	r = np.abs(fermions.x1 - fermions.x2)
	r = np.where(r < 0.0001, 0.0001, r)
	V = k/ r
	return V

H = Hamiltonian(particles = TwoFermions(), 
				potential = coulomb_interaction, # change this to harmonic_interaction to check what happens!
				spatial_ndim = 1, N = 200, extent = 10*Å)


eigenstates = H.solve(max_states = 90)
print(eigenstates.energies)

visualization = init_visualization(eigenstates)

#visualization.slider_plot()
visualization.animate()