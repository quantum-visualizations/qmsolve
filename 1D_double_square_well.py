import numpy as np
from qmsolve import Hamiltonian,  SingleParticle, init_visualization


def double_well(particle):
    k = 50
    N = len(particle.x)
    W, U, L = 20, 11, 9
    V_lst = [1.0 if (i < N//10 or i > 9*N//10 or 
                     (i > L*N//W  and i < U*N//W))
             else 0.0 for i in range(N)]
    return 2.0*k*np.array(V_lst)



H = Hamiltonian(particles = SingleParticle(), 
				potential = double_well, 
				spatial_ndim = 1, N = 100, extent = 20)


eigenstates = H.solve(max_states = 10)

vis = init_visualization(eigenstates)
coeffs = np.zeros([10], np.complex128)
coeffs[0], coeffs[1] = 1.0, 1.0j
vis.superpositions(states=coeffs, dt=300.0)