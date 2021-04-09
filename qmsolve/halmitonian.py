import numpy as np
from scipy.sparse.linalg import eigsh
from scipy.sparse import diags
import time


class Halmitonian:
    def __init__(self,particles, potential, N, extent, spatial_ndim):
        """
        N: number of grid points
        extent: spacial extent, measured in angstroms
        """

        
        self.N = N
        self.extent = extent
        self.dx = extent/N
        self.particle_system = particles
        self.spatial_ndim = spatial_ndim 
        self.ndim = 0 # total number of observables 

        self.particle_system.get_observables(self)
        self.T = self.particle_system.get_kinetic_matrix(self)

        self.potential = potential
        self.V = self.get_potential_matrix()

        
    def get_potential_matrix(self):

        V = self.potential(self.particle_system)
        V = V.reshape(self.N**self.ndim)
        V = diags([V], [0])
        return V

    def solve(self, max_states):

        H = self.T + self.V
        print ("Computing...")
        t0 = time.time()
        eigenvalues, eigenvectors = eigsh(H, k = max_states, which='SA')

        """
        the result of this method depends of the particle system. For example if the systems are two fermions, this method
        makes the eigenstates antisymetric
        """
        self.energies, self.eigenstates = self.particle_system.get_energies_and_eigenstates(self, max_states, eigenvalues, eigenvectors)

        print ("Took", time.time() - t0)

        return self.energies, self.eigenstates
