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


L = 20
N = 200

# well for the xy direction
def xy_well(particle):
    N = len(particle.x)
    V = np.zeros([N, N])
    height = 150.0
    V[:, 64*N//128: 65*N//128] = height
    V[64*N//128: 65*N//128, :] = height
    V[0: 4*N//32, :] = height
    V[:, 0: 4*N//32] = height
    V[N - 4*N//32: N, :] = height
    V[:, N - 4*N//32: N] = height
    return V


# well for the z direction
def z_well(particle):
    N = len(particle.x)
    V = np.zeros([N])
    height = 150.0
    V[64*N//128: 65*N//128] = height
    V[0: 4*N//32] = height
    V[N - 4*N//32: N] = height
    return V

M_XY = 9
M_Z = 3
H_xy = Hamiltonian(particles = SingleParticle(), 
				   potential = xy_well, 
				   spatial_ndim = 2, N = N, extent = L)
e_xy, v_xy = H_xy.solve(max_states = M_XY)
# dynamic_visualize(e_xy, v_xy)
H_z = Hamiltonian(particles = SingleParticle(),
				  potential = z_well,
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