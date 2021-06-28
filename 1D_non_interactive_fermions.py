import numpy as np
from qmsolve import Hamiltonian, TwoFermions, init_visualization, Å


H = Hamiltonian(particles = TwoFermions(), 
				potential = None, # If None is specified, the particles are just limited by the limits of the grid which are infinite barriers
				spatial_ndim = 1, N = 100, extent = 10 * Å)


eigenstates = H.solve(max_states = 90)
print(eigenstates.energies)

visualization = init_visualization(eigenstates)

#visualization.plot_eigenstate(6)
#visualization.slider_plot()
visualization.animate()