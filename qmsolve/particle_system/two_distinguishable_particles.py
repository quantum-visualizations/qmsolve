import numpy as np
from scipy.sparse import diags
from scipy.sparse import kron
from scipy.sparse import eye
from .two_particles import TwoParticles
from ..util.constants import *
from .. import Eigenstates

class TwoDistinguishableParticles(TwoParticles):

    def __init__(self, m1 = m_e, m2 = m_e, spin = None):
        """
        N: number of grid points
        extent: spacial extent, measured in angstroms
        """
        self.m1 = m1
        self.m2 = m2

        self.spin = spin

    def get_kinetic_matrix(self, H):

        I = eye(H.N)
        T_ =  diags([-2., 1., 1.], [0,-1, 1] , shape=(H.N, H.N))*-k/(H.dx**2)

        if H.spatial_ndim ==1:
            T =  (kron(T_/self.m1,I) + kron(I,T_/self.m2))
        elif H.spatial_ndim ==2:
            T =  (kron(T_/self.m1,I,I,I) + kron(I,T_/self.m1,I,I) + kron(I,I,T_/self.m2,I) + kron(I,I,I,T_/self.m2))

        return T



    def get_energies_and_eigenstates(self, H, max_states, eigenvalues, eigenvectors):

        energies = eigenvalues
        eigenstates = eigenvectors.T.reshape(( max_states, *[H.N]*H.ndim) )

        # Finish the normalization of the eigenstates
        eigenstates = eigenstates/np.sqrt(H.dx**H.ndim)
        return  energies/eV, eigenstates