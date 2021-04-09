import numpy as np
from scipy.sparse import diags
from scipy.sparse import kron
from scipy.sparse import eye
from .particle_system import Particle_system
from ..util.constants import *


class Two_fermions(Particle_system):
    def __init__(self, m = m_e, spin = None):
        """
        N: number of grid points
        extent: spacial extent, measured in angstroms
        """
        self.m = m
        self.spin = spin


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



    def get_kinetic_matrix(self, H):

        I = eye(H.N)
        T_ =  diags([-2., 1., 1.], [0,-1, 1] , shape=(H.N, H.N))*-k/(self.m*H.dx**2)

        if H.spatial_ndim ==1:
            T =  (kron(T_,I) + kron(I,T_))
        elif H.spatial_ndim ==2:
            T =  (kron(T_,I,I,I) + kron(I,T_,I,I) + kron(I,I,T_,I) + kron(I,I,I,T_))

        return T


    #Not tested yet. Work in progress:
    def get_energies_and_eigenstates(self, H, max_states, eigenvalues, eigenvectors):

        eigenvectors = eigenvectors = eigenvectors.T.reshape(( max_states, *[H.N]*H.ndim) )


        energies = []
        antisymmetric_eigenvectors = []

        #antisymmetrize eigenvectors:
        for i in range(max_states):
            for j in range(i+1, max_states):
                if i != j:
                    energies +=  [eigenvalues[i] + eigenvalues[j]]
                    antisymmetric_eigenvectors +=  [eigenvectors[i]*eigenvectors[j] - eigenvectors[j]*eigenvectors[i]]

        #sort energies
        def take_eigenvalues(lst):
            return lst[0]
        lst  = list(zip(eigenvalues,eigenvectors))
        lst.sort(key=take_eigenvalues)
        energies,eigenstates = zip(*lst)
        eigenstates = np.array(eigenstates)


        # Normalize the eigenstates

        return energies, eigenstates