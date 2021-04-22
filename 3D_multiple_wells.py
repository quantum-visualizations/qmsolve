"""
3D well, where the potential takes the separable form
V(x, y, z) = V_{xy}(x, y) + V_{z}(z)
"""

import numpy as np
from qmsolve import Hamiltonian, SingleParticle
try:
    from qmsolve import visualize3D, animate3D, dynamic_visualize
except ImportError as e:
    print("You do not have Mayavi.")
    raise e


L = 10
N = 100

# well for the xy direction
def four_gaussian_wells(particle):
	𝜇 = 2.0
	σ = 0.6
	V = 300*(3-np.exp((-(particle.x)**2 -(particle.y-𝜇)**2 ) / (2*σ**2)) / (np.sqrt(2*np.pi)* σ)
	-np.exp((-(particle.x-𝜇)**2 -(particle.y)**2 ) / (2*σ**2)) / (np.sqrt(2*np.pi)* σ)
	-np.exp((-(particle.x+𝜇)**2 -(particle.y)**2 ) / (2*σ**2)) / (np.sqrt(2*np.pi)* σ)
	-np.exp((-(particle.x)**2 -(particle.y+𝜇)**2 ) / (2*σ**2)) / (np.sqrt(2*np.pi)* σ))
	return V

# well for the z direction
def double_well(particle):
	return 50.0*(particle.x**2 - 1.0**2)**2

M_XY = 9
M_Z = 3
H_xy = Hamiltonian(particles = SingleParticle(), 
				   potential = four_gaussian_wells, 
				   spatial_ndim = 2, N = N, extent = L)
e_xy, v_xy = H_xy.solve(max_states = M_XY)
# dynamic_visualize(e_xy, v_xy)
H_z = Hamiltonian(particles = SingleParticle(),
				  potential = double_well,
				  spatial_ndim = 1, N = N, extent=L)
e_z, v_z = H_z.solve(max_states = M_Z)
# dynamic_visualize(e_z, v_z)


states = np.transpose(np.multiply.outer(v_xy, v_z), (0, 3, 1, 2, 4))
states = states.reshape([M_XY*M_Z, N, N, N])
energies = np.add.outer(e_xy, e_z).reshape([M_XY*M_Z])
sort_array = np.argsort(energies)
energies = energies[sort_array]
states = states[sort_array]

visualize3D(energies, states, 0, plot_type='contour')
animate3D(energies, states)