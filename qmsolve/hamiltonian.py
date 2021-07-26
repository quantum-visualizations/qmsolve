import numpy as np
from scipy.sparse import diags
import time


class Hamiltonian:
    def __init__(self, particles, potential, N, extent, spatial_ndim, potential_type = "grid", E_min=0):
        """Class to represent the Hamiltonian operator.

        Parameters
        ----------
        particles : ParticleSystem
            The particle system
        potential : callable
            The potential
        N : int
            number of grid points
        extent : array-like
            spacial extent, measured in angstroms
        spatial_ndim : [type]
            Initial guess for the energy of the ground state. It's only used if potential_type = "matrix" is used
        potential_type : str, optional
            The type of the potential, by default "grid"
        E_min : int, optional
            The minimum energy, by default 0
        """

        self.N = N
        self.extent = extent
        self.dx = extent / N
        self.particle_system = particles
        self.spatial_ndim = spatial_ndim
        self.ndim = 0  # total number of observables

        self.T = self.particle_system.get_kinetic_matrix(self)

        self.potential = potential
        self.potential_type = potential_type
        self.E_min = E_min

        if potential_type == "grid":
            self.particle_system.get_observables(self)

        elif potential_type == "matrix":
            self.particle_system.build_matrix_operators(self)

        self.V = self.get_potential_matrix()

    def get_potential_matrix(self):
        """Get the potential matrix.

        Returns
        -------
        array-like
            The potential matrix
        """

        if self.potential_type == "grid":

            if self.potential == None:
                self.E_min = 0.
                V = 0.
                return V
            else: 
                V = self.potential(self.particle_system)
                self.E_min = np.amin(V)
                V = V.reshape(self.N ** self.ndim)
                V = diags([V], [0])
                return V

        elif self.potential_type == "matrix":
             return self.potential(self.particle_system)



    def solve(self, max_states: int, method: str = 'eigsh', N0 = 30, maxiter = 30, verbose = False):
        """Solve the Hamiltonian for its energies and eigenstates.

        Parameters
        ----------
        max_states : int
            The maximum number of states to solve.
        method : str
            Which methods to use (currently only "eigsh", "lobpcg", and "lobpcg-cupy" are available), by default "eigsh".
        N0 : int
            TODO
        maxiter : int
            The number of iterations to perform. by default False
        verbose : bool
            TODO
        """
        implemented_solvers = ('eigsh', 'lobpcg', 'lobpcg-cupy')

        H = self.T + self.V
        print("Computing...")

        t0 = time.time()

        if method == 'eigsh':
            from scipy.sparse.linalg import eigsh

            # Note: uses shift-invert trick for stability finding low-lying states
            # Ref: https://docs.scipy.org/doc/scipy/reference/tutorial/arpack.html#shift-invert-mode

            eigenvalues, eigenvectors = eigsh(H, k=max_states, which='LM', sigma=min(0, self.E_min))



        elif method == 'lobpcg':
            from scipy.sparse.linalg import eigsh, lobpcg, LinearOperator
            from scipy.sparse import diags
            if self.spatial_ndim != 3:
                raise NotImplementedError(
                    "lobpcg is only implemented for a 3D single particle")

            from qmsolve import SingleParticle
            #First, we compute eighs eigenvectors with a grid of size N0, 
            H_eigsh = Hamiltonian(particles = SingleParticle(), 
                                  potential = self.potential, 
                                  spatial_ndim = 3, N = N0, extent = self.extent, potential_type = self.potential_type, E_min = self.E_min)

            eigenvalues_eigsh, eigenvectors_eigsh = eigsh(H_eigsh.V + H_eigsh.T, k=max_states, which='LM', sigma=min(0, self.E_min))

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
                    "lobpcg is only implemented for a 3D single particle")

            from qmsolve import SingleParticle
            #First, we compute eighs eigenvectors with a grid of size N0, 
            H_eigsh = Hamiltonian(particles = SingleParticle(), 
                                  potential = self.potential, 
                                  spatial_ndim = 3, N = N0, extent = self.extent, potential_type = self.potential_type, E_min = self.E_min)

            eigenvalues_eigsh, eigenvectors_eigsh = eigsh(H_eigsh.V + H_eigsh.T, k=max_states, which='LM', sigma=min(0, self.E_min))

            eigenvectors_eigsh = eigenvectors_eigsh.reshape(  *[N0]*3 , max_states)

            if verbose == True:
                print("Initial eigsh computation completed")

            if self.potential_type == "grid":
                #Now, we interpolate them to a grid of size N and then use it as an initial guess to the lobpcg solver.
                from scipy.interpolate import interpn
                new_xx, new_yy, new_zz, states = np.mgrid[ -1:1:self.N*1j, -1:1:self.N*1j, -1:1:self.N*1j, -1:1:max_states*1j]
                eigenvectors_eigsh_interpolated = interpn((np.linspace(-1,1,N0), np.linspace(-1,1,N0), np.linspace(-1,1,N0), np.linspace(-1,1,max_states)), 
                                                          eigenvectors_eigsh, 
                                                          np.array([new_xx, new_yy, new_zz, states]).T).T

            elif self.potential_type == "matrix":
                raise NotImplementedError(
                 "lobpcg-cupy solver has not been implemented to work with complex numbers. Use lobpcg instead")



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
                str(method) + " solver has not been implemented. Use one of " + str(implemented_solvers))

        """the result of this method depends of the particle system. For example if the systems are two fermions, 
        this method makes the eigenstates antisymmetric """
        self.eigenstates = self.particle_system.get_eigenstates(self, max_states, eigenvalues, eigenvectors)

        # When using complex numbers in the potential energies aren't necessarily sorted
        if self.potential_type == "matrix":
            sort_array = np.argsort(self.eigenstates.energies)
            self.eigenstates.energies = self.eigenstates.energies[sort_array]
            self.eigenstates.array = self.eigenstates.array[sort_array]


        print("Took", time.time() - t0)
        return self.eigenstates
