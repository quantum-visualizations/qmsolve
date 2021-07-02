import numpy as np
from qmsolve import Hamiltonian, SingleParticle, init_visualization,Å


#interaction potential
def four_gaussian_wells(particle):
	𝜇 = 1.6*Å
	σ = 0.5*Å
	V = 22*(4-np.exp((-(particle.x)**2 -(particle.y-𝜇)**2 ) / (2*σ**2))
	-np.exp((-(particle.x-𝜇)**2 -(particle.y)**2 ) / (2*σ**2))
	-np.exp((-(particle.x+𝜇)**2 -(particle.y)**2 ) / (2*σ**2))
	-np.exp((-(particle.x)**2 -(particle.y+𝜇)**2 ) / (2*σ**2)))
	return V



H = Hamiltonian(particles = SingleParticle(), 
				potential = four_gaussian_wells, 
				spatial_ndim = 2, N = 200, extent = 8*Å)


eigenstates = H.solve(max_states = 60)

print(eigenstates.energies)
visualization = init_visualization(eigenstates)
visualization.animate()
coeffs = np.zeros([30], np.complex128)
coeffs[6] = 1.0
coeffs[7] = 1.0j
coeffs[28] = (1.0 - 1.0j)/np.sqrt(2.0)
coeffs[29] = (1.0 + 1.0j)/np.sqrt(2.0)
visualization.superpositions(coeffs, dt=0.0002, 
							 # hide_controls=True, 
							 # save_animation=True, frames=7*30
							 )


