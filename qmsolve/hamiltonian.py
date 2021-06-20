import numpy as np
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


    def solve(self, max_states: int, method: str = 'eigsh', N0 = 30, maxiter = 30, verbose = False):
        """
        Diagonalize the hamiltonian and retrieve the lowest-energy eigenstates
        Args:
            max_states: the number of states to retreive
            method: the solver method. Currently, 'eigsh' and 'lobpcg' are implemented. Note: 'lobpcg' is potentially
            much faster than 'eigsh' but can fail catastrophically for some systems. Use 'lobpcg' with care.
            N0: grid divisions for the initial eigsh computations. Later the eigenstates will be scales
        Returns:

        """
        implemented_solvers = ('eigsh', 'lobpcg', 'lobpcg-cupy')

        H = self.T + self.V
        print("Computing...")

        t0 = time.time()

        if method == 'eigsh':
            from scipy.sparse.linalg import eigsh

            # Note: uses shift-invert trick for stability finding low-lying states
            # Ref: https://docs.scipy.org/doc/scipy/reference/tutorial/arpack.html#shift-invert-mode

            eigenvalues, eigenvectors = eigsh(H, k=max_states, which='LM', sigma=min(0, self.Vmin))



        elif method == 'lobpcg':
            from scipy.sparse.linalg import eigsh, lobpcg, LinearOperator
            from scipy.sparse import diags
            if self.spatial_ndim != 3:
                raise NotImplementedError(
                    f"lobpcg is only implemented for a 3D single particle")

            from qmsolve import SingleParticle
            #First, we compute eighs eigenvectors with a grid of size N0, 
            H_eigsh = Hamiltonian(particles = SingleParticle(), 
                                  potential = self.potential, 
                                  spatial_ndim = 3, N = N0, extent = self.extent)

            eigenvalues_eigsh, eigenvectors_eigsh = eigsh(H_eigsh.V + H_eigsh.T, k=max_states, which='LM', sigma=min(0, self.Vmin))

            eigenvectors_eigsh = eigenvectors_eigsh.reshape(  *[N0]*3 , max_states)

            if verbose == True:
                print("Initial eigsh computation completed")


            #Now, we interpolate them to a grid of size N and then use it as an initial guess to the lobpcg solver.
            from scipy.interpolate import interpn
            new_xx, new_yy, new_zz, states = np.mgrid[ -1:1:self.N*1j, -1:1:self.N*1j, -1:1:self.N*1j, -1:1:max_states*1j]
            eigenvectors_eigsh_interpolated = interpn((np.linspace(-1,1,N0), np.linspace(-1,1,N0), np.linspace(-1,1,N0), np.linspace(-1,1,max_states)), 
                                                      eigenvectors_eigsh, 
                                                      np.array([new_xx, new_yy, new_zz, states]).T).T
            if verbose == True:
                print("Interpolation completed")

            # preconditioning matrix should approximate the inverse of the hamiltonian
            # we naively construct this by taking the inverse of diagonal elements
            # and setting all others to zero. This is called the Jacobi or diagonal preconditioner.
            A = diags([1 / H.diagonal()], [0])
            precond = lambda x: A @ x
            M = LinearOperator(H.shape, matvec=precond, matmat=precond)
            # guess for eigenvectors is computed from random numbers
            # TODO: a better guess would be to use the eigenstates of a reference system like a square well
            eigenvectors_guess = eigenvectors_eigsh_interpolated.reshape(  self.N**self.ndim , max_states)

            sol = lobpcg(H, eigenvectors_guess, largest=False, M=M, tol=1e-15, maxiter = maxiter)
            eigenvalues, eigenvectors = sol[0], sol[1]

            if verbose == True:
                print("lobpcg computation completed")

        elif method == 'lobpcg-cupy':
            from scipy.sparse.linalg import eigsh

            if self.spatial_ndim != 3:
                raise NotImplementedError(
                    f"lobpcg is only implemented for a 3D single particle")

            from qmsolve import SingleParticle
            #First, we compute eighs eigenvectors with a grid of size N0, 
            H_eigsh = Hamiltonian(particles = SingleParticle(), 
                                  potential = self.potential, 
                                  spatial_ndim = 3, N = N0, extent = self.extent)

            eigenvalues_eigsh, eigenvectors_eigsh = eigsh(H_eigsh.V + H_eigsh.T, k=max_states, which='LM', sigma=min(0, self.Vmin))

            eigenvectors_eigsh = eigenvectors_eigsh.reshape(  *[N0]*3 , max_states)

            if verbose == True:
                print("Initial eigsh computation completed")

            #Now, we interpolate them to a grid of size N and then use it as an initial guess to the lobpcg solver.
            from scipy.interpolate import interpn
            new_xx, new_yy, new_zz, states = np.mgrid[ -1:1:self.N*1j, -1:1:self.N*1j, -1:1:self.N*1j, -1:1:max_states*1j]
            eigenvectors_eigsh_interpolated = interpn((np.linspace(-1,1,N0), np.linspace(-1,1,N0), np.linspace(-1,1,N0), np.linspace(-1,1,max_states)), 
                                                      eigenvectors_eigsh, 
                                                      np.array([new_xx, new_yy, new_zz, states]).T).T
            if verbose == True:
                print("Interpolation completed")
            eigenvectors_guess = eigenvectors_eigsh_interpolated.reshape(  self.N**self.ndim , max_states)

            from cupyx.scipy.sparse.linalg import lobpcg, LinearOperator
            from cupyx.scipy.sparse import diags
            from cupyx.scipy.sparse.csr import csr_matrix
            # preconditioning matrix should approximate the inverse of the hamiltonian
            # we naively construct this by taking the inverse of diagonal elements
            # and setting all others to zero. This is called the Jacobi or diagonal preconditioner.
            H = csr_matrix(H)
            A = diags([1 / H.diagonal()], [0]).tocsc()
            precond = lambda x: A @ x
            M = LinearOperator(H.shape, matvec=precond, matmat=precond)

            # guess for eigenvectors is computed from random numbers
            # TODO: a better guess would be to use the eigenstates of a reference system like a square well
            import cupy as cp
            sol = lobpcg(H, cp.array(eigenvectors_guess), largest=False, M=M, tol=1e-15, maxiter = maxiter)
            eigenvalues, eigenvectors = sol[0].get(), sol[1].get()


        else:
            raise NotImplementedError(
                f"{method} solver has not been implemented. Use one of {implemented_solvers}")

        """the result of this method depends of the particle system. For example if the systems are two fermions, 
        this method makes the eigenstates antisymmetric """
        self.eigenstates = self.particle_system.get_eigenstates(self, max_states, eigenvalues, eigenvectors)

        print("Took", time.time() - t0)
        return self.eigenstates
