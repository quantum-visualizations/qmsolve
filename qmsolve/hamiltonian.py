import numpy as np
from scipy.sparse import diags
import time


class Hamiltonian:
    def __init__(self, particles, potential, N, extent, spatial_ndim, potential_type = "grid", E_min=0):
        """
        N: number of grid points
        extent: spacial extent, measured in bohr radius (length atomic unit)
        E_min: Initial guess for the energy of the ground state measured in hartrees (energy atomic unit). It's only used if potential_type = "matrix" is used
        """

        self.N = N
        self.extent = extent
        self.dx = extent / N
        self.particle_system = particles
        self.particle_system.H = self
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

        if self.potential_type == "grid":

            if self.potential == None:
                self.E_min = 0.
                V = 0.
                return V
            else: 
                V = self.potential(self.particle_system)
                self.Vgrid = V
                self.E_min = np.amin(V)
                V = V.reshape(self.N ** self.ndim)
                V = diags([V], [0])
                return V

        elif self.potential_type == "matrix":
             V =  self.potential(self.particle_system)
             self.Vgrid = np.real((V).diagonal().reshape(*([self.N] *self.ndim )) + self.E_min)

             # Note: Vgrid when potential_type == "matrix" is only used for visualization. 
             # It represents the potential without the effect of momentum terms
             return V





    def solve(self, max_states: int, method: str = 'eigsh', verbose = False, lobpcg_args = {'N0': 30, 'preconditioner' : 'jacobi', 'maxiter' : 30}):
        """
        Diagonalize the hamiltonian and retrieve the lowest-energy eigenstates
        Args:
            max_states: the number of states to retreive
            method: the solver method. Currently, 'eigsh' and 'lobpcg' are implemented.
            lobpcg_args:
                N0: grid divisions for the initial eigsh computations to be used as an initial guess in lobpcg.
                preconditioner: lobpcg preconditioner. 'pyamg' convergence is faster but requires having installed pyamg and may fail for some hamiltonians.
                                Default preconditioner is 'jacobi'.
                maxiter: maximum number of iterations. 
        Returns:
            eigenstates
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
            implemented_lobpcg_preconditioners = ('jacobi', 'pyamg')

            if self.spatial_ndim != 3:
                raise NotImplementedError(
                    f"lobpcg is only implemented for a 3D single particle")

            from qmsolve import SingleParticle
            #First, we compute eighs eigenvectors with a grid of size N0, 
            H_eigsh = Hamiltonian(particles = SingleParticle(), 
                                  potential = self.potential, 
                                  spatial_ndim = 3, N = lobpcg_args['N0'], extent = self.extent, potential_type = self.potential_type, E_min = self.E_min)

            eigenvalues_eigsh, eigenvectors_eigsh = eigsh(H_eigsh.V + H_eigsh.T, k=max_states, which='LM', sigma=min(0, self.E_min))

            eigenvectors_eigsh = eigenvectors_eigsh.reshape(  *[lobpcg_args['N0']]*3 , max_states)

            if verbose == True:
                print("Initial eigsh computation completed")


            #Now, we interpolate them to a grid of size N and then use it as an initial guess to the lobpcg solver.
            from scipy.interpolate import interpn
            new_xx, new_yy, new_zz, states = np.mgrid[ -1:1:self.N*1j, -1:1:self.N*1j, -1:1:self.N*1j, -1:1:max_states*1j]
            eigenvectors_eigsh_interpolated = interpn((np.linspace(-1,1,lobpcg_args['N0']), np.linspace(-1,1,lobpcg_args['N0']), np.linspace(-1,1,lobpcg_args['N0']), np.linspace(-1,1,max_states)), 
                                                      eigenvectors_eigsh, 
                                                      np.array([new_xx, new_yy, new_zz, states]).T).T
            eigenvectors_guess = eigenvectors_eigsh_interpolated.reshape(  self.N**self.ndim , max_states)

            if verbose == True:
                print("Interpolation completed")

            
            if lobpcg_args['preconditioner'] == 'jacobi':
                # preconditioning matrix should approximate the inverse of the hamiltonian
                # we naively construct this by taking the inverse of diagonal elements
                # and setting all others to zero. This is called the Jacobi or diagonal preconditioner.
                A = diags([1 / H.diagonal()], [0])
                precond = lambda x: A @ x
                M = LinearOperator(H.shape, matvec=precond, matmat=precond)

            elif lobpcg_args['preconditioner'] == 'pyamg':
                # to install pyamg run 'pip install pyamg'
                from pyamg import smoothed_aggregation_solver
                ml = smoothed_aggregation_solver(H)
                M = ml.aspreconditioner()
            else:
                raise NotImplementedError(
                    f"{lobpcg_args['preconditioner']} preconditioner has not been implemented. Use one of {implemented_lobpcg_preconditioners}")


            sol = lobpcg(H, eigenvectors_guess, largest=False, M=M, tol=1e-15, maxiter = lobpcg_args['maxiter'])
            eigenvalues, eigenvectors = sol[0], sol[1]

            if verbose == True:
                print("lobpcg computation completed")

        elif method == 'lobpcg-cupy':
            from scipy.sparse.linalg import eigsh
            implemented_lobpcg_preconditioners = ('jacobi')

            if self.spatial_ndim != 3:
                raise NotImplementedError(
                    f"lobpcg is only implemented for a 3D single particle")

            from qmsolve import SingleParticle
            #First, we compute eighs eigenvectors with a grid of size N0, 
            H_eigsh = Hamiltonian(particles = SingleParticle(), 
                                  potential = self.potential, 
                                  spatial_ndim = 3, N = lobpcg_args['N0'], extent = self.extent, potential_type = self.potential_type, E_min = self.E_min)

            eigenvalues_eigsh, eigenvectors_eigsh = eigsh(H_eigsh.V + H_eigsh.T, k=max_states, which='LM', sigma=min(0, self.E_min))

            eigenvectors_eigsh = eigenvectors_eigsh.reshape(  *[lobpcg_args['N0']]*3 , max_states)

            if verbose == True:
                print("Initial eigsh computation completed")

            if self.potential_type == "grid":
                #Now, we interpolate them to a grid of size N and then use it as an initial guess to the lobpcg solver.
                from scipy.interpolate import interpn
                new_xx, new_yy, new_zz, states = np.mgrid[ -1:1:self.N*1j, -1:1:self.N*1j, -1:1:self.N*1j, -1:1:max_states*1j]
                eigenvectors_eigsh_interpolated = interpn((np.linspace(-1,1,lobpcg_args['N0']), np.linspace(-1,1,lobpcg_args['N0']), np.linspace(-1,1,lobpcg_args['N0']), np.linspace(-1,1,max_states)), 
                                                          eigenvectors_eigsh, 
                                                          np.array([new_xx, new_yy, new_zz, states]).T).T

            elif self.potential_type == "matrix":
                raise NotImplementedError(
                f"lobpcg-cupy solver has not been implemented to work with complex numbers. Use lobpcg instead")



            if verbose == True:
                print("Interpolation completed")
            eigenvectors_guess = eigenvectors_eigsh_interpolated.reshape(self.N**self.ndim , max_states)

            from cupyx.scipy.sparse.linalg import lobpcg, LinearOperator
            from cupyx.scipy.sparse import diags
            from cupyx.scipy.sparse.csr import csr_matrix
            H = csr_matrix(H)

            if lobpcg_args['preconditioner'] == 'jacobi':
                # preconditioning matrix should approximate the inverse of the hamiltonian
                # we naively construct this by taking the inverse of diagonal elements
                # and setting all others to zero. This is called the Jacobi or diagonal preconditioner.
                A = diags([1 / H.diagonal()], [0]).tocsc()
                precond = lambda x: A @ x
                M = LinearOperator(H.shape, matvec=precond, matmat=precond)
            else:
                raise NotImplementedError(
                    f"{lobpcg_args['preconditioner']} preconditioner has not been implemented. Use one of {implemented_lobpcg_preconditioners}")

            import cupy as cp
            sol = lobpcg(H, cp.array(eigenvectors_guess), largest=False, M=M, tol=1e-15, maxiter = lobpcg_args['maxiter'])
            eigenvalues, eigenvectors = sol[0].get(), sol[1].get()


        else:
            raise NotImplementedError(
                f"{method} solver has not been implemented. Use one of {implemented_solvers}")

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
