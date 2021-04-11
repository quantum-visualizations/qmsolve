import numpy as np
from qmsolve import Halmitonian, animate, dynamic_visualize, SingleParticle


#interaction potential
def four_gaussian_wells(particle):
	𝜇 = 2
	σ = 0.6
	V = 200*(3-np.exp((-(particle.x)**2 -(particle.y-𝜇)**2 ) / (2*σ**2)) / (np.sqrt(2*np.pi)* σ)
	-np.exp((-(particle.x-𝜇)**2 -(particle.y)**2 ) / (2*σ**2)) / (np.sqrt(2*np.pi)* σ)
	-np.exp((-(particle.x+𝜇)**2 -(particle.y)**2 ) / (2*σ**2)) / (np.sqrt(2*np.pi)* σ)
	-np.exp((-(particle.x)**2 -(particle.y+𝜇)**2 ) / (2*σ**2)) / (np.sqrt(2*np.pi)* σ))
	return V



H = Halmitonian(particles = SingleParticle(), 
				potential = four_gaussian_wells, 
				spatial_ndim = 2, N = 100, extent = 10)


energies, eigenstates = H.solve(max_states = 30)


animate(energies, eigenstates)
dynamic_visualize(energies, eigenstates)