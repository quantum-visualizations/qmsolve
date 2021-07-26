import numpy as np
from qmsolve import visualization
from qmsolve import Hamiltonian, SingleParticle, Eigenstates,Å


#interaction potential
def harmonic_oscillator(particle):

	kx, ky, kz = [0.020]*3 # measured in eV / (Å**2)

	return 0.500 * kx * particle.x**2 +\
		   0.51 * ky * particle.y**2 +\
		   0.52 * kz * particle.z**2


H = Hamiltonian(particles=SingleParticle(),
                  N=90, potential=harmonic_oscillator, 
                  extent=12*Å, spatial_ndim=3)


eigenstates = H.solve( max_states = 32, method ='lobpcg')
print(eigenstates.energies)


visualization = visualization.init_visualization(eigenstates)
visualization.plot_type = 'contour'

#visualization.plot_eigenstate(0)

#visualization.plot_eigenstate(26)
visualization.animate()
