import numpy as np
from qmsolve import Hamiltonian, SingleParticle
from qmsolve import animate, dynamic_visualize, visualize
from qmsolve import visualize_superpositions


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


energies, eigenstates = H.solve(max_states = 10)

print(energies)
visualize_superpositions(energies, eigenstates, 10, dt=0.01)