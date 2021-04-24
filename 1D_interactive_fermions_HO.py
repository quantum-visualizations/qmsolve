import numpy as np
from qmsolve import Hamiltonian, TwoFermions, init_visualization


#interaction potential
def harmonic_oscillator_plus_coulomb_interaction(fermions):

	k = 100 # measured in eV / (Å**2)

	V_harmonic = 0.5*k*fermions.x1**2 + 0.5*k*fermions.x2**2 

	k = 300. # measured in eV * Å**2
	r = np.abs(fermions.x1 - fermions.x2)
	r = np.where(r < 0.0001, 0.0001, r)
	V_interaction = k/ r

	return V_harmonic + V_interaction




H = Hamiltonian(particles = TwoFermions(), 
				potential = harmonic_oscillator_plus_coulomb_interaction, 
				spatial_ndim = 1, N = 200, extent = 15)


eigenstates = H.solve(max_states = 32)
print(eigenstates.energies)

visualization = init_visualization(eigenstates)
visualization.animate(max_states = 32, xlim = [-4,4])