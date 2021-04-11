import numpy as np
from qmsolve import Halmitonian, TwoFermions, dynamic_visualize, animate


def harmonic_interaction(fermions):

	k = 50 # measured in eV / (Å**2)

	l0 = 3 # measured in Å
	V = 0.5*k*(fermions.x1 - fermions.x2 - l0) **2
	return V



def coulomb_interaction(fermions):

	k = 140. # measured in eV * Å**2
	r = (fermions.x1 - fermions.x2)
	r = np.where(r < 0.0001, 0.0001, r)
	V = k/ r**2
	return V

H = Halmitonian(particles = TwoFermions(), 
				potential = harmonic_interaction, # change this to coulomb_interaction to check what happens!
				spatial_ndim = 1, N = 100, extent = 10)


energies, eigenstates = H.solve(max_states = 90)
print("Energies:",energies)
#dynamic_visualize(energies, eigenstates)
animate(energies, eigenstates)
