import numpy as np
from scipy.sparse import diags
from scipy.sparse import kron
from scipy.sparse import eye
from .particle_system import ParticleSystem
from ..util.constants import *


class SingleParticle(ParticleSystem):
    def __init__(self, m = m_e, spin = None):
        """
        N: number of grid points
        extent: spacial extent, measured in angstroms
        """
        self.m = m
        self.spin = spin


    def get_observables(self, H):

        if H.spatial_ndim ==1:
            self.x = np.linspace(-H.extent/2, H.extent/2, H.N)
            H.ndim = 1

        elif H.spatial_ndim ==2:
            x = np.linspace(-H.extent/2, H.extent/2, H.N)
            y = np.linspace(-H.extent/2, H.extent/2, H.N)
            self.x, self.y = np.meshgrid(x,y)
            H.ndim = 2


        elif H.spatial_ndim ==3:
            x = np.linspace(-H.extent/2, H.extent/2, H.N)
            y = np.linspace(-H.extent/2, H.extent/2, H.N)
            z = np.linspace(-H.extent/2, H.extent/2, H.N)
            self.x, self.y, self.z  = np.meshgrid(x,y,z)
            H.ndim = 3

    def get_kinetic_matrix(self, H):

        I = eye(H.N)
        T_ =  diags([-2., 1., 1.], [0,-1, 1] , shape=(H.N, H.N))*-k/(self.m*H.dx**2)
        if H.spatial_ndim ==1:
            T = T_

        elif H.spatial_ndim ==2:
            T =  (kron(T_,I) + kron(I,T_))

        elif H.spatial_ndim ==3:
            T =  (kron(T_,I,I) + kron(I,T_,I) + kron(I,I,T_))

        return T

    def get_energies_and_eigenstates(self, H, max_states, eigenvalues, eigenvectors):

        energies = eigenvalues
        eigenstates = eigenvectors.T.reshape(( max_states, *[H.N]*H.ndim) )

        # Finish the normalization of the eigenstates
        eigenstates = eigenstates/np.sqrt(H.dx**H.ndim)
        return  energies, eigenstates