import numpy as np
from qmsolve import Halmitonian, dynamic_visualize


H = Halmitonian(N = 80, extent = 10)
H.add_particle(spatial_ndim = 2)


L = 10
N = 80
x = np.linspace(-L/2,L/2,N)
y = np.linspace(-L/2,L/2,N)
xx, yy = np.meshgrid(x, y)
𝜇 = 2
σ = 0.5
V = 200*(3-np.exp((-(xx)**2 -(yy-𝜇)**2 ) / (2*σ**2)) / (np.sqrt(2*np.pi)* σ)
-np.exp((-(xx-𝜇)**2 -(yy)**2 ) / (2*σ**2)) / (np.sqrt(2*np.pi)* σ)
-np.exp((-(xx+𝜇)**2 -(yy)**2 ) / (2*σ**2)) / (np.sqrt(2*np.pi)* σ)
-np.exp((-(xx)**2 -(yy+𝜇)**2 ) / (2*σ**2)) / (np.sqrt(2*np.pi)* σ))


H.add_potential(V)
energies, eigenstates = H.solve(30)


#visualize(energies, eigenstates, 23) #static version, only plots a eigenstate
dynamic_visualize(energies, eigenstates)