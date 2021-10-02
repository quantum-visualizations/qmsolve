import numpy as np
from .method import Method
import time
from ..util.constants import *
from ..particle_system import SingleParticle, TwoParticles
from scipy import sparse
from scipy.sparse import linalg
import progressbar

"""
Crank-Nicolson method for the Schrödinger equation: https://imsc.uni-graz.at/haasegu/Lectures/HPC-II/SS17/presentation1_Schroedinger-Equation_HPC2-seminar.pdf
Prototype and original implementation: https://gist.github.com/marl0ny/23947165652ccad73e55b01241afbe77 

Jacobi iteration can be used to solve for the system of equations that occurs
when using the Crank-Nicolson method. This is following an
article found here https://arxiv.org/pdf/1409.8340 by Sadovskyy et al.,
where the Crank-Nicolson method with Jacobi iteration is 
used to solve the Ginzburg-Landau equations, which is similar to 
the Schrödinger equation but contains nonlinear terms and
couplings with the four vector potential.
"""



def jacobi(inv_diag: sparse.dia_matrix, lower_upper: sparse.dia.dia_matrix,
           b: np.ndarray, min_iter: int = 10,  max_iter: int = 20, TOL = 0.001):
    """
    Given the inverse diagonals of a matrix A and its lower and upper parts
    without its diagonal, find the solution x for the system Ax = b.
    Reference for Jacobi iteration:
    https://en.wikipedia.org/wiki/Jacobi_method
    """
    x = b.copy()
    for i in range(min_iter):
        x = inv_diag @ (b - lower_upper @ x)
        
    for i in range(max_iter - min_iter):
        x_ = inv_diag @ (b - lower_upper @ x)
        rel_err = np.mean(np.abs(x - x_))
        x = x_
        if rel_err < TOL:
            break
    return x



class CrankNicolson(Method):
    def __init__(self, simulation):

        self.simulation = simulation
        self.H = simulation.H

        if self.H.potential_type == "matrix":
            self.H.particle_system.get_observables(self.H)

        self.simulation.Vmin = np.amin(self.H.Vgrid)
        self.simulation.Vmax = np.amax(self.H.Vgrid)

    def run(self, initial_wavefunction, total_time, dt, store_steps = 1):

        self.simulation.store_steps = store_steps
        dt_store = total_time/store_steps
        self.simulation.total_time = total_time

        Nt = int(np.round(total_time / dt))
        Nt_per_store_step = int(np.round(dt_store / dt))
        self.simulation.Nt_per_store_step = Nt_per_store_step

        #time/dt and dt_store/dt must be integers. Otherwise dt is rounded to match that the Nt_per_store_stepdivisions are integers
        self.simulation.dt = dt_store/Nt_per_store_step

        if isinstance(self.simulation.H.particle_system ,SingleParticle):
            Ψ = np.zeros((store_steps + 1, self.H.N **self.H.ndim), dtype = np.complex128)
            I = sparse.identity(self.H.N **self.H.ndim)
            Ψ[0] = np.array(initial_wavefunction(self.H.particle_system)).reshape( self.H.N **self.H.ndim)

        elif isinstance(self.simulation.H.particle_system,TwoParticles):
            Ψ = np.zeros((store_steps + 1, self.H.N ** 2), dtype = np.complex128)
            I = sparse.identity(self.H.N ** 2)
            Ψ[0] = np.array(initial_wavefunction(self.H.particle_system)).reshape(self.H.N**2 )




        m = self.H.particle_system.m


        BETA = 0.5j*self.simulation.dt/hbar

        H_matrix = self.H.T + self.H.V
        A = I + BETA*H_matrix
        B = I - BETA*H_matrix
        #We are going to solve the equation A*Ψ_{i+1} = B*Ψ_{i} for Ψ_{i+1}

        D = sparse.diags(A.diagonal(0), (0))
        INV_D = sparse.diags(1.0/A.diagonal(0), (0))
        L_PLUS_U = A - D


        t0 = time.time()
        bar = progressbar.ProgressBar()
        for i in bar(range(store_steps)):
            tmp = np.copy(Ψ[i])
            for j in range(Nt_per_store_step):
                B_dot_Ψ = B @ tmp
                tmp = linalg.gcrotmk(A, B_dot_Ψ)[0]
                #tmp = jacobi(INV_D, L_PLUS_U, B_dot_Ψ, min_iter=10, max_iter = 50, TOL = 0.00001)
            Ψ[i+1] = tmp
        print("Took", time.time() - t0)

        if isinstance(self.simulation.H.particle_system ,SingleParticle):
            self.simulation.Ψ = Ψ.reshape(store_steps + 1, *([self.H.N] *self.H.ndim ))

        elif isinstance(self.simulation.H.particle_system,TwoParticles):
            self.simulation.Ψ = Ψ.reshape(store_steps + 1, *([self.H.N] *2 ))

        self.simulation.Ψmax = np.amax(np.abs(Ψ))









class CrankNicolsonCupy(Method):
    def __init__(self, simulation):


        self.simulation = simulation
        self.H = simulation.H

        if self.H.potential_type == "matrix":
            self.H.particle_system.get_observables(self.H)

        self.simulation.Vmin = np.amin(self.H.Vgrid)
        self.simulation.Vmax = np.amax(self.H.Vgrid)

    def run(self, initial_wavefunction, total_time, dt, store_steps = 1):
        import cupy as cp 
        from cupyx.scipy import sparse

        self.simulation.store_steps = store_steps
        dt_store = total_time/store_steps
        self.simulation.total_time = total_time

        Nt = int(np.round(total_time / dt))
        Nt_per_store_step = int(np.round(dt_store / dt))
        self.simulation.Nt_per_store_step = Nt_per_store_step

        #time/dt and dt_store/dt must be integers. Otherwise dt is rounded to match that the Nt_per_store_stepdivisions are integers
        self.simulation.dt = dt_store/Nt_per_store_step

        if isinstance(self.simulation.H.particle_system ,SingleParticle):
            Ψ = cp.zeros((store_steps + 1, self.H.N **self.H.ndim), dtype = cp.complex128)
            I = sparse.identity(self.H.N **self.H.ndim)
            Ψ[0] = cp.array(initial_wavefunction(self.H.particle_system)).reshape( self.H.N **self.H.ndim)

        elif isinstance(self.simulation.H.particle_system,TwoParticles):
            Ψ = cp.zeros((store_steps + 1, self.H.N ** 2), dtype = cp.complex128)
            I = sparse.identity(self.H.N ** 2)
            Ψ[0] = cp.array(initial_wavefunction(self.H.particle_system)).reshape(self.H.N**2 )




        m = self.H.particle_system.m


        BETA = 0.5j*self.simulation.dt/hbar

        H_matrix = sparse.csr.csr_matrix(self.H.T + self.H.V)
        A = I + BETA*H_matrix
        B = I - BETA*H_matrix
        #We are going to solve the equation A*Ψ_{i+1} = B*Ψ_{i} for Ψ_{i+1}

        D = sparse.diags(A.diagonal(0), (0))
        INV_D = sparse.diags(1.0/A.diagonal(0), (0))
        L_PLUS_U = A - D


        def jacobi_cupy(inv_diag: sparse.dia_matrix, lower_upper: sparse.dia.dia_matrix,
                   b: cp.ndarray, min_iter: int = 10,  max_iter: int = 20, TOL = 0.001):
            """
            Given the inverse diagonals of a matrix A and its lower and upper parts
            without its diagonal, find the solution x for the system Ax = b.
            Reference for Jacobi iteration:
            https://en.wikipedia.org/wiki/Jacobi_method
            """
            x = b.copy()
            for i in range(min_iter):
                x = inv_diag @ (b - lower_upper @ x)
                
            for i in range(max_iter - min_iter):
                x_ = inv_diag @ (b - lower_upper @ x)
                rel_err = cp.mean(cp.abs(x - x_))
                x = x_
                if rel_err < TOL:
                    break
            return x


        bar = progressbar.ProgressBar()
        t0 = time.time()
        for i in bar(range(store_steps)):
            tmp = cp.copy(Ψ[i])
            for j in range(Nt_per_store_step):
                B_dot_Ψ = B @ tmp
                tmp = jacobi_cupy(INV_D, L_PLUS_U, B_dot_Ψ, min_iter=10, max_iter = 50, TOL = 0.00001)
            Ψ[i+1] = tmp
        print("Took", time.time() - t0)

        Ψ = Ψ.get()
        if isinstance(self.simulation.H.particle_system ,SingleParticle):
            self.simulation.Ψ = Ψ.reshape(store_steps + 1, *([self.H.N] *self.H.ndim ))

        elif isinstance(self.simulation.H.particle_system,TwoParticles):
            self.simulation.Ψ = Ψ.reshape(store_steps + 1, *([self.H.N] *2 ))

        self.simulation.Ψmax = np.amax(np.abs(Ψ))


