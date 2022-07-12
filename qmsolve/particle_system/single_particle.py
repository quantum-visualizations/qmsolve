import numpy as np
from scipy.sparse import diags
from scipy.sparse import kron
from scipy.sparse import eye
from .particle_system import ParticleSystem
from ..util.constants import *
from .. import Eigenstates

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
            self.x, self.y, self.z  = np.mgrid[ -H.extent/2: H.extent/2:H.N*1j, -H.extent/2: H.extent/2:H.N*1j, -H.extent/2: H.extent/2:H.N*1j]
            H.ndim = 3

    def compute_momentum_space(self, H):
        """
        Used for split step method
        """

        if H.spatial_ndim == 1:

            self.px = np.fft.fftshift(np.fft.fftfreq(H.N, d = H.dx)) * hbar  * 2*np.pi
            self.py = np.fft.fftshift(np.fft.fftfreq(H.N, d = H.dx)) * hbar  * 2*np.pi
            self.p2 = self.px**2


        elif H.spatial_ndim ==2:

            px = np.fft.fftshift(np.fft.fftfreq(H.N, d = H.dx)) * hbar  * 2*np.pi
            py = np.fft.fftshift(np.fft.fftfreq(H.N, d = H.dx)) * hbar  * 2*np.pi
            px, py = np.meshgrid(px, py)


            self.p2 = (px**2 + py**2)

        elif self.H.spatial_ndim == 3:
            raise NotImplementedError(
                f"split-step isn't implemented for a 3D single particle")



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
            

            self.I = kron(I,kron(I,I))

    def get_kinetic_matrix(self, H):

        I = eye(H.N)
        T_ =  diags([-2., 1., 1.], [0,-1, 1] , shape=(H.N, H.N))*-k/(self.m*H.dx**2)
        if H.spatial_ndim ==1:
            T = T_

        elif H.spatial_ndim ==2:
            T =  kron(T_,I) + kron(I,T_)

        elif H.spatial_ndim ==3:
            T =  kron(T_, kron(I, I)) + kron(I, kron(T_, I)) + kron(I, kron(I, T_))

        return T

    def get_eigenstates(self, H, max_states, eigenvalues, eigenvectors):

        energies = eigenvalues
        eigenstates_array = np.moveaxis(eigenvectors.reshape(  *[H.N]*H.ndim , max_states), -1, 0)

        # Finish the normalization of the eigenstates
        eigenstates_array = eigenstates_array/np.sqrt(H.dx**H.ndim)

        if H.spatial_ndim == 1:
            type = "SingleParticle1D"
        elif H.spatial_ndim == 2:
            type = "SingleParticle2D"
        elif H.spatial_ndim == 3:
            type = "SingleParticle3D"

        eigenstates = Eigenstates(energies/eV, eigenstates_array, H.extent, H.N, type)
        return eigenstates