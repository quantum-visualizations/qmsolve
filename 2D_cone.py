import numpy as np
from qmsolve import Hamiltonian, animate, dynamic_visualize
from qmsolve import SingleParticle, visualize_superpositions



def cone(particle):
    return 5.0*np.sqrt(particle.x**2 + particle.y**2)



H = Hamiltonian(particles = SingleParticle(), 
				potential = cone, 
				spatial_ndim = 2, N = 100, extent = 15)


energies, eigenstates = H.solve(max_states = 11)

print(energies)
visualize_superpositions(energies, eigenstates, 10, dt=0.01)
dynamic_visualize(energies, eigenstates)
