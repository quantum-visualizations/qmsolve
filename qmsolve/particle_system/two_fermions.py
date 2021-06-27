import numpy as np
from scipy.sparse import diags
from scipy.sparse import kron
from scipy.sparse import eye
from .two_particles import TwoParticles
from ..util.constants import *
from .. import Eigenstates


class TwoFermions(TwoParticles):


    def get_eigenstates(self, H, max_states, eigenvalues, eigenvectors):

        eigenvectors  = eigenvectors.T.reshape(( max_states, *[H.N]*H.ndim) )

        # Normalize the eigenvectors
        eigenvectors = eigenvectors/np.sqrt(H.dx**H.ndim)
        

        energies = []
        eigenstates_array = []

        #antisymmetrize eigenvectors: This is made by applying (𝜓(r1 , s1, r2 , s2) - 𝜓(r2 , s2, r1 , s1))/sqrt(2) to each state.
        for i in range(max_states):
            eigenstate_tmp = (eigenvectors[i] - eigenvectors[i].swapaxes(0,1))/np.sqrt(2)

            norm = np.sum(eigenstate_tmp*eigenstate_tmp)*H.dx**H.ndim 

            TOL = 0.02
            
            # check if is eigenstate_tmp is a normalizable eigenstate. (norm shouldn't be zero)
            if norm > TOL : 
                # for some reason when the eigenstate is degenerated it isn't normalized 
                #print("norm",norm)
                eigenstate_tmp = eigenstate_tmp/np.sqrt(norm)

                 
                if eigenstates_array != []: #check if it's the first eigenstate
                    inner_product = np.sum(eigenstates_array[-1]* eigenstate_tmp)*H.dx**H.ndim
                    #print("inner_product",inner_product)
                else:
                    inner_product = 0


                if np.abs(inner_product) < TOL: # check if is eigenstate_tmp is repeated. (inner_product should be zero)

                    eigenstates_array +=  [eigenstate_tmp]
                    energies +=  [eigenvalues[i]]

        if H.spatial_ndim == 1:
            type = "TwoIdenticalParticles1D"
        elif H.spatial_ndim == 2:
            type = "TwoIdenticalParticles2D"

        eigenstates = Eigenstates(np.array(energies)/eV, eigenstates_array, H.extent, H.N, type)
        return eigenstates