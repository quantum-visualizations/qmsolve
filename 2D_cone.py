import numpy as np
from qmsolve import Hamiltonian, SingleParticle, init_visualization



def cone(particle):
    return 5.0*np.sqrt(particle.x**2 + particle.y**2)



H = Hamiltonian(particles = SingleParticle(), 
				potential = cone, 
				spatial_ndim = 2, N = 100, extent = 15)


eigenstates = H.solve(max_states = 30)

visualization = init_visualization(eigenstates)
#visualization.plot_eigenstate(6)
coeffs = np.zeros([30], np.complex128)
coeffs[6] = 1.0
coeffs[7] = 1.0j
coeffs[28] = (1.0 - 1.0j)/np.sqrt(2.0)
coeffs[29] = (1.0 + 1.0j)/np.sqrt(2.0)
visualization.superpositions(coeffs, dt=0.002, 
							 # hide_controls=True, 
							 # save_animation=True, frames=7*30
							 )
# visualization.animate()

