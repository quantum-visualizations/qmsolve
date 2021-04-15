import numpy as np
from scipy.sparse import diags
from scipy.sparse import kron
from scipy.sparse import eye
from .particle_system import ParticleSystem
from ..util.constants import *


class TwoIdenticalParticles(ParticleSystem):


    def __init__(self, m = m_e, spin = None, p=1):
        """
        N: number of grid points
        extent: spacial extent, measured in angstroms
        """
        self._p = p
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


    #Not fully tested yet.
    def get_energies_and_eigenstates(self, H, max_states, eigenvalues, eigenvectors):

        eigenvectors  = eigenvectors.T.reshape(( max_states, *[H.N]*H.ndim) )

        # Normalize the eigenvectors
        eigenvectors = eigenvectors/np.sqrt(H.dx**H.ndim)
        

        energies = []
        eigenstates = []

        #antisymmetrize eigenvectors: This is made by applying (𝜓(r1 , s1, r2 , s2) - 𝜓(r2 , s2, r1 , s1))/sqrt(2) to each state.
        for i in range(max_states):
            eigenstate_tmp = (eigenvectors[i] + self._p*eigenvectors[i].swapaxes(0,1))/np.sqrt(2)

            norm = np.sum(eigenstate_tmp*eigenstate_tmp)*H.dx**H.ndim 

            TOL = 0.02
            
            # check if is eigenstate_tmp is a normalizable eigenstate. (norm shouldn't be zero)
            if norm > TOL : 
                # for some reason when the eigenstate is degenerated it isn't normalized 
                #print("norm",norm)
                eigenstate_tmp = eigenstate_tmp/np.sqrt(norm)

                 
                if eigenstates != []: #check if it's the first eigenstate
                    inner_product = np.sum(eigenstates[-1]* eigenstate_tmp)*H.dx**H.ndim
                    #print("inner_product",inner_product)
                else:
                    inner_product = 0


                if np.abs(inner_product) < TOL: # check if is eigenstate_tmp is repeated. (inner_product should be zero)

                    eigenstates +=  [eigenstate_tmp]
                    energies +=  [eigenvalues[i]]

        return energies, eigenstates
        