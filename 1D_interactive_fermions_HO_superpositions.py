import numpy as np
from qmsolve import Hamiltonian, TwoFermions, init_visualization


#interaction potential
def harmonic_oscillator_plus_coulomb_interaction(fermions):

	k = 100 # measured in eV / (Å**2)

	V_harmonic = 0.5*k*fermions.x1**2 + 0.5*k*fermions.x2**2 

	k = 300. # measured in eV * Å
	r = np.abs(fermions.x1 - fermions.x2)
	r = np.where(r < 0.0001, 0.0001, r)
	V_interaction = k/ r

	return V_harmonic + V_interaction




H = Hamiltonian(particles = TwoFermions(), 
				potential = harmonic_oscillator_plus_coulomb_interaction, 
				spatial_ndim = 1, N = 100, extent = 15)


eigenstates = H.solve(max_states = 200)

visualization = init_visualization(eigenstates)


#This example visualizes the time dependent Schrödinger of two interactive fermions with a a superposition of the eigenstates
x1, x2 = np.meshgrid(np.linspace(-7.5, 7.5, eigenstates.array[0].shape[0]),
                np.linspace(-7.5, 7.5, eigenstates.array[0].shape[1]))

mu0 = 2.0
𝜓0 = np.exp(-(x1+mu0)**2/(2*0.375**2))*np.exp(-(x2-mu0)**2/(2*0.375**2)) - np.exp(-(x1-mu0)**2/(2*0.375**2))*np.exp(-(x2+mu0)**2/(2*0.375**2))

#compute the inner product of the initial state 𝜓0(x1,x2) with the eigenstates 𝜓_n(x1,x2):  
#coeffs = <𝜓_n(x1,x2)|𝜓0(x1,x2)>
coeffs = np.tensordot(eigenstates.array, 𝜓0, axes=([1, 2], [0,1]))*1.0j
#visualize a superposition of the eigenstates
visualization.superpositions(coeffs,xlim=[-3.5, 3.5], hide_controls = True)
