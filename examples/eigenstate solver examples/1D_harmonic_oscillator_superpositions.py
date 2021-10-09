import numpy as np
from qmsolve import Hamiltonian,  SingleParticle, init_visualization, Å,eV


#=========================================================================================================#
#We define the Hamiltonian of a single particle confined in an harmonic oscillator potential. 
#Then, we compute its eigenstates.
#=========================================================================================================#


#interaction potential
def harmonic_oscillator(particle):

	k = 100 * eV / Å**2
	return 0.5 * k * particle.x**2


#define the Hamiltonian
H = Hamiltonian(particles = SingleParticle(), 
				potential = harmonic_oscillator, 
				spatial_ndim = 1, N = 512, extent = 20*Å)

#Diagonalize the Hamiltonian and compute the eigenstates
eigenstates = H.solve(max_states = 30)

#=========================================================================================================#
#The next lines are used for visualizing a superposition of eigenstates, where 𝜓0 is a gaussian wave packet.
#=========================================================================================================#


x = np.linspace(-1.0*Å, 1.0*Å, len(eigenstates.array[0]))
𝜓0 = np.exp(-(x-0.16*Å)**2/(2*(0.05*Å)**2))

#compute the inner product of the initial state 𝜓0(x) with the eigenstates 𝜓_n(x). (coeffs = <𝜓_n|𝜓0>)
coeffs = np.dot(eigenstates.array, 𝜓0)*1.0j

#visualize a superposition of the eigenstates
visualization = init_visualization(eigenstates)
visualization.superpositions(coeffs[0:15], xlim=[-3.5*Å, 3.5*Å])
