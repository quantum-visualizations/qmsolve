import numpy as np
from scipy.sparse import diags
from scipy.sparse import kron
from scipy.sparse import eye
from .particle_system import ParticleSystem
from ..util.constants import *
from abc import abstractmethod 

class TwoParticles(ParticleSystem):

    @abstractmethod
    def __init__(self, m = m_e, spin = None):
        """
        N: number of grid points
        extent: spacial extent, measured in angstroms
        """
        self.m = m
        self.spin = spin

    @abstractmethod
    def get_observables(self, H):

        if H.spatial_ndim ==1:
            x1 = np.linspace(-H.extent/2, H.extent/2, H.N)
            x2 = np.linspace(-H.extent/2, H.extent/2, H.N)
            self.x1, self.x2 = np.meshgrid(x1,x2)
            H.ndim = 2

        elif H.spatial_ndim ==2:
            x1 = np.linspace(-H.extent/2, H.extent/2, H.N)
            y1 = np.linspace(-H.extent/2, H.extent/2, H.N)
            x2 = np.linspace(-H.extent/2, H.extent/2, H.N)
            y2 = np.linspace(-H.extent/2, H.extent/2, H.N)
            H.ndim = 4

            self.x1, self.y1, self.x2, self.y2 = np.meshgrid(x1,y1,x2,y2)


    @abstractmethod
    def get_kinetic_matrix(self, H):

        I = eye(H.N)
        T_ =  diags([-2., 1., 1.], [0,-1, 1] , shape=(H.N, H.N))*-k/(self.m*H.dx**2)

        if H.spatial_ndim ==1:
            T =  (kron(T_,I) + kron(I,T_))
        elif H.spatial_ndim ==2:
            T =  (kron(T_,I,I,I) + kron(I,T_,I,I) + kron(I,I,T_,I) + kron(I,I,I,T_))

        return T



        