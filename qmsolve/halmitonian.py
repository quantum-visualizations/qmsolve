import numpy as np
from scipy.sparse.linalg import eigsh
from scipy.sparse import diags
from scipy.sparse import kron
from scipy.sparse import eye

import matplotlib.pyplot as plt
import time

k = 3.8099821161548593 # hbar**2 / (2*m_e) /(Å**2) / eV
m_e = 1



class Halmitonian:
    def __init__(self, N, extent):
        """
        N: number of grid points
        extent: spacial extent, measured in angstroms
        """

        
        self.N = N
        self.extent = extent
        self.dx = extent/N
        self.ndim = 0

    def add_particle(self, spatial_ndim, m = m_e, spin = None):
        # prototype version. If we are going to solve for several particles, we should have a separate Particle class.

        T_ =  diags([-2., 1., 1.], [0,-1, 1] , shape=(self.N, self.N))*-k/(m*self.dx**2)
        I = eye(self.N)

        if spatial_ndim ==1:
            T = T_
            self.ndim = 1

        elif spatial_ndim ==2:
            T =  (kron(T_,I) + kron(I,T_))
            self.ndim = 2


        elif spatial_ndim ==3:
            T =  (kron(T_,I,I) + kron(I,T_,I) + kron(I,I,T_))
            self.ndim = 3

        self.T = T
        
    def add_potential(self, V):

        self.V = V.reshape(self.N**self.ndim)
        self.V = diags([self.V], [0])

        self.H = self.V + self.T

    def solve(self, max_states):

        H = self.T + self.V
        print ("Computing...")
        t0 = time.time()
        self.energies, self.eigenstates = eigsh(H, k = max_states, which='SA')
        self.eigenstates = self.eigenstates.T.reshape(( max_states, *[self.N]*self.ndim) )
        print ("Took", time.time() - t0)

        # Finish the normalization of the eigenstates
        self.eigenstates = self.eigenstates/(self.dx**self.ndim)

        return self.energies, self.eigenstates
