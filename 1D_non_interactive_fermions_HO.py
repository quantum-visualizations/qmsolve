import numpy as np
from qmsolve import Hamiltonian, TwoFermions, init_visualization


#interaction potential
def harmonic_oscillator(fermions):

	k = 100 # measured in eV / (Å**2)

	V = 0.5*k*fermions.x1**2 + 0.5*k*fermions.x2**2 
	return V



H = Hamiltonian(particles = TwoFermions(), 
				potential = harmonic_oscillator, 
				spatial_ndim = 1, N = 200, extent = 15)


eigenstates = H.solve(max_states = 80)
print(eigenstates.energies)

visualization = init_visualization(eigenstates)
visualization.animate(max_states = 31)