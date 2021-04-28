import numpy as np
from qmsolve import Hamiltonian, SingleParticle, init_visualization


#interaction potential
def two_gaussian_wells(particle):
	𝜇 = 0.7
	σ = 0.5
	V = 700*(2+
	-np.exp((-(particle.x-𝜇)**2 -(particle.y)**2 ) / (2*σ**2))
	-np.exp((-(particle.x+𝜇)**2 -(particle.y)**2 ) / (2*σ**2)))
	return V



H = Hamiltonian(particles = SingleParticle(), 
				potential = two_gaussian_wells, 
				spatial_ndim = 2, N = 100, extent = 8)


eigenstates = H.solve(max_states = 40)

print(eigenstates.energies)
visualization = init_visualization(eigenstates)

visualization.animate()

#visualize a superposition of eigenstates
coeffs = np.zeros([10], np.complex128)
coeffs[0] = 1.0
coeffs[1] = -1.0
visualization.superpositions(coeffs, dt=0.03, xlim=[-3.0, 3.0],
							 # save_animation=True, frames=60
							 )
