import numpy as np
from scipy.sparse import diags
from scipy.sparse import kron
from scipy.sparse import eye
from .particle_system import ParticleSystem
from ..util.constants import *
from .. import Eigenstates

class SingleParticle(ParticleSystem):
    def __init__(self, m = m_e, spin = None, symmetry = None, azimuthal_number_m = 0):
        """
        N: number of grid points
        extent: spacial extent, measured in angstroms
        """
        self.m = m
        self.spin = spin
        self.symmetry = symmetry
        self.azimuthal_number_m = azimuthal_number_m

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
            if self.symmetry == "azimuthal":
                r = np.linspace(0., H.extent/2, H.N+1)[1:]
                z = np.linspace(-H.extent/2, H.extent/2, H.N)
                #φ = np.linspace(0, 2*np.pi, H.N)
                self.r, self.z = np.meshgrid(r, z)
                H.ndim = 2

            else:
                self.x, self.y, self.z  = np.mgrid[ -H.extent/2: H.extent/2:H.N*1j, -H.extent/2: H.extent/2:H.N*1j, -H.extent/2: H.extent/2:H.N*1j]
                H.ndim = 3



    def build_matrix_operators(self, H):

        if H.spatial_ndim == 1:

            H.ndim = 1
            
            self.x = np.linspace(-H.extent/2, H.extent/2, H.N)
            diff_x =  diags([-1., 0., 1.], [-1, 0, 1] , shape=(H.N, H.N))*1/(2*H.dx)
            self.px = - hbar *1j * diff_x
            
            self.I = eye(H.N)


        elif H.spatial_ndim == 2:

            H.ndim = 2

            x = diags([np.linspace(-H.extent/2, H.extent/2, H.N)], [0])
            y = diags([np.linspace(-H.extent/2, H.extent/2, H.N)], [0])
            I = eye(H.N)

            self.x = kron(I,x)
            self.y = kron(y,I)

            diff_x = diags([-1., 0., 1.], [-1, 0, 1] , shape=(H.N, H.N))*1/(2*H.dx)
            diff_y = diags([-1., 0., 1.], [-1, 0, 1] , shape=(H.N, H.N))*1/(2*H.dx)

            self.px = kron(I, - hbar *1j * diff_y)
            self.py = kron(- hbar *1j * diff_x, I)
            
            self.I = kron(I,I)

        elif H.spatial_ndim == 3:

            H.ndim = 3
            xx, yy, zz  = np.mgrid[ -H.extent/2: H.extent/2:H.N*1j, -H.extent/2: H.extent/2:H.N*1j, -H.extent/2: H.extent/2:H.N*1j]
            r = np.sqrt(xx**2 + yy**2 + zz**2)

            TOL = 0.000001
            r_inv = np.where(r < TOL, 1/TOL, 1./r)

            self.r = diags([ r.reshape(H.N ** H.ndim) ], [0])   
            self.r_inv = diags([ r_inv.reshape(H.N ** H.ndim) ], [0])   

            x = diags([np.linspace(-H.extent/2, H.extent/2, H.N)], [0])
            y = diags([np.linspace(-H.extent/2, H.extent/2, H.N)], [0])
            z = diags([np.linspace(-H.extent/2, H.extent/2, H.N)], [0])
            I = eye(H.N)

            self.x = kron(x, kron(I, I))
            self.y = kron(I, kron(y, I))
            self.z = kron(I, kron(I, z))

            diff_x = diags([-1., 0., 1.], [-1, 0, 1] , shape=(H.N, H.N))*1/(2*H.dx)
            diff_y = diags([-1., 0., 1.], [-1, 0, 1] , shape=(H.N, H.N))*1/(2*H.dx)
            diff_z = diags([-1., 0., 1.], [-1, 0, 1] , shape=(H.N, H.N))*1/(2*H.dx)

            self.px = kron(- hbar *1j * diff_x, kron(I,  I))
            self.py = kron(I , kron(- hbar *1j * diff_x, I))
            self.pz = kron(I, kron(I , - hbar *1j * diff_x))

    def get_kinetic_matrix(self, H):

        I = eye(H.N)
        T_ =  diags([-2., 1., 1.], [0,-1, 1] , shape=(H.N, H.N))*-k/(self.m*H.dx**2)
        if H.spatial_ndim ==1:
            T = T_

        elif H.spatial_ndim ==2:
            T =  kron(T_,I) + kron(I,T_)

        elif H.spatial_ndim ==3:
            if self.symmetry == "azimuthal":
                Tz = T_

                r = np.linspace(0., H.extent/2, H.N+1)[1:]
                self.dr = r[1]-r[0]
                H.dr = self.dr
                
                diag1 = np.append(1, np.ones(H.N-1))
                diag2 = np.append(1, -2*np.ones(H.N-1))
                diag3 = np.append(-2,  np.ones(H.N-1))
                diag4 = np.append(1,  np.zeros(H.N-1))
                Dr2 = diags([diag1, diag2,  diag3, diag4]  ,  [-1, 0, 1,2], shape=(H.N, H.N))*1/(H.dr**2)

                diag1 = -np.ones(H.N)
                diag2 = np.append(-2, np.zeros(H.N-1))
                diag3 = np.append(2, np.ones(H.N-1))
                Dr = diags([diag1, diag2,  diag3]  ,  [-1, 0, 1], shape=(H.N, H.N))*1/(2*H.dr)

                
                sqrt_r = diags([np.sqrt(r)], [0])  # r factor
                H.r = r
                inv_r = diags([1./r], [0])  # 1/r factor
                inv_r2 = diags([1./(r**2)], [0]) # 1 / (r**2) factor

                T =  kron(Tz, I) + kron( I , -k/self.m *(Dr2 + inv_r * Dr  - inv_r2*self.azimuthal_number_m**2) )
            else:

                T =  kron(T_, kron(I, I)) + kron(I, kron(T_, I)) + kron(I, kron(I, T_))

        return T

    def get_eigenstates(self, H, max_states, eigenvalues, eigenvectors):

        energies = eigenvalues
        eigenstates_array = np.moveaxis(eigenvectors.reshape(  *[H.N]*H.ndim , max_states), -1, 0)

        # Finish the normalization of the eigenstates

        if self.symmetry == "azimuthal":

            r_fact = ( np.outer(np.ones(max_states), self.r)).reshape(max_states,H.N,H.N)
            norm = np.sum(eigenstates_array*eigenstates_array * r_fact)*H.dr*2*np.pi
            eigenstates_array /np.sqrt(norm)
        else:
            eigenstates_array = eigenstates_array/np.sqrt(H.dx**H.ndim)

        if H.spatial_ndim == 1:
            type = "SingleParticle1D"
        elif H.spatial_ndim == 2:
            type = "SingleParticle2D"
        elif H.spatial_ndim == 3:
            type = "SingleParticle3D"

        eigenstates = Eigenstates(energies/eV, eigenstates_array, H.extent, H.N, type)
        return eigenstates