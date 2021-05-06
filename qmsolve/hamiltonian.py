import numpy as np
from scipy.sparse.linalg import eigsh, lobpcg, LinearOperator
from scipy.sparse import diags
import time


class Hamiltonian:
    def __init__(self, particles, potential, N, extent, spatial_ndim):
        """
        N: number of grid points
        extent: spacial extent, measured in angstroms
        """

        self.N = N
        self.extent = extent
        self.dx = extent / N
        self.particle_system = particles
        self.spatial_ndim = spatial_ndim
        self.ndim = 0  # total number of observables

        self.particle_system.get_observables(self)
        self.T = self.particle_system.get_kinetic_matrix(self)

        self.potential = potential
        self.V = self.get_potential_matrix()

    def get_potential_matrix(self):
        if self.potential == None:
            self.Vmin = 0.
            V = 0.
            return V
        else: 
            V = self.potential(self.particle_system)
            self.Vmin = np.amin(V)
            V = V.reshape(self.N ** self.ndim)
            V = diags([V], [0])
            return V


    def solve(self, max_states: int, method: str = 'eigsh'):
        """
        Diagonalize the hamiltonian and retrieve the lowest-energy eigenstates
        Args:
            max_states: the number of states to retreive
            method: the solver method. Currently, 'eigsh' and 'lobpcg' are implemented. Note: 'lobpcg' is potentially
            much faster than 'eigsh' but can fail catastrophically for some systems. Use 'lobpcg' with care.

        Returns:

        """
        implemented_solvers = ('eigsh', 'lobpcg')

        H = self.T + self.V
        print("Computing...")
        t0 = time.time()

        if method == 'eigsh':
            # Note: uses shift-invert trick for stability finding low-lying states
            # Ref: https://docs.scipy.org/doc/scipy/reference/tutorial/arpack.html#shift-invert-mode

            eigenvalues, eigenvectors = eigsh(H, k=max_states, which='LM', sigma=min(0, self.Vmin))
        elif method == 'lobpcg':
            # preconditioning matrix should approximate the inverse of the hamiltonian
            # we naively construct this by taking the inverse of diagonal elements
            # and setting all others to zero. This is called the Jacobi or diagonal preconditioner.
            A = diags([1 / H.diagonal()], [0]).tocsc()
            precond = lambda x: A @ x
            M = LinearOperator(H.shape, matvec=precond, matmat=precond)

            # guess for eigenvectors is computed from random numbers
            # TODO: a better guess would be to use the eigenstates of a reference system like a square well
            X_approx = np.random.rand(H.shape[0], max_states)

            sol = lobpcg(H, X_approx, largest=False, M=M, tol=1e-15)
            eigenvalues, eigenvectors = sol[0], sol[1]
        else:
            raise NotImplementedError(
                f"{method} solver has not been implemented. Use one of {implemented_solvers}")

        """the result of this method depends of the particle system. For example if the systems are two fermions, 
        this method makes the eigenstates antisymmetric """
        self.eigenstates = self.particle_system.get_eigenstates(self, max_states, eigenvalues, eigenvectors)

        print("Took", time.time() - t0)
        return self.eigenstates
