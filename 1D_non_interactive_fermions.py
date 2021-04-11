import numpy as np
from qmsolve import Halmitonian, TwoFermions, dynamic_visualize, animate


#interaction potential
def harmonic_oscillator(fermions):

	k = 100 # measured in eV / (Å**2)

	V = 0.5*k*fermions.x1**2 + 0.5*k*fermions.x2**2 
	return V



H = Halmitonian(particles = TwoFermions(), 
				potential = harmonic_oscillator, 
				spatial_ndim = 1, N = 300, extent = 15)


energies, eigenstates = H.solve(max_states = 30)
print("Energies:",energies)
animate(energies, eigenstates)
#dynamic_visualize(energies, eigenstates)
