import numpy as np
from qmsolve import visualization
from qmsolve import Hamiltonian, SingleParticle, Eigenstates


spatial_ndim = 1
L = 10
N = 100
separable_potential = [lambda particle: 20.0*particle.x**2,
                       lambda particle: 20.0*particle.x**2,
                       lambda particle: 20.0*particle.x**2]

energies = []
states = []
M = 3 # Number of 1D eigenstates to solve.

for V_i in separable_potential:
    H_i = Hamiltonian(particles=SingleParticle(),
                      N=N, potential=V_i, 
                      extent=L, spatial_ndim=1)
    v = H_i.solve(M)
    energies.append(v.energies)
    states.append(v.array)

states = np.multiply.outer(np.multiply.outer(states[2], states[1]), states[0])
states = np.transpose(states, (0, 2, 4, 1, 3, 5))
energies = np.add.outer(np.add.outer(energies[2], energies[1]), energies[0])
energies = energies.flatten()
states = states.reshape([M*M*M, N, N, N])
sort_array = np.argsort(energies)
# How to rearrange an array with argsort output:
# https://stackoverflow.com/a/64638032
energies = energies[sort_array]
states = states[sort_array]


eigenstates = Eigenstates(energies, states, H_i, 'SingleParticle3D')
v = visualization.init_visualization(eigenstates)
v.plot_type = 'contour'
v.superpositions(np.array([0.0]*10 + [((1.0 + 1.0j)/np.sqrt(2.0))**i for i in range(10)]), 
                 dt=0.005)
v.plot_type = 'volume'
v.animate()
