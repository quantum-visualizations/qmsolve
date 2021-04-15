"""
Visualization of the numerically computed
eigenstates of the 3D harmonic oscillator (WIP).
Separation of variables is used for faster computation.
"""
import numpy as np
from qmsolve import Hamiltonian, SingleParticle
try:
    from qmsolve import visualize3D, animate3D
except ImportError as e:
    print("You do not have Mayavi.")
    raise e


spatial_ndim = 1
L = 20
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
    e_and_v = H_i.solve(M)
    energies.append(e_and_v[0])
    states.append(e_and_v[1])

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

visualize3D(energies, states, 1)
# Set delay from its default value of 500 to something smaller for a faster 
# animation.
animate3D(energies, states)
