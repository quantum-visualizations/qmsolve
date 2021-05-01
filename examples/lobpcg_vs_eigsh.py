import numpy as np
from qmsolve import Hamiltonian, SingleParticle, init_visualization


#interaction potential
def harmonic_oscillator(particle):

	kx, ky, kz = [2]*3 # measured in eV / (Å**2)
	return 0.5 * kx * particle.x**2 +\
		   0.5 * ky * particle.y**2 +\
		   0.5 * kz * particle.z**2


if __name__=='__main__':

	# simulation params
	grid_size=30
	extent=10
	num_states=10

	# set up the system
	H = Hamiltonian(particles = SingleParticle(),
					potential = harmonic_oscillator,
					spatial_ndim = 3, N = grid_size, extent = extent)


	# solve with eigsh
	print('-'*30)
	print("Solving with scipy.sparse.linalg.eigsh")
	eig_eigsh = H.solve(max_states = num_states).energies
	print("Eigenvalues (eigsh) :", eig_eigsh)

	# solve with lobpcg
	print('-'*30)
	print("Solving with scipy.sparse.linalg.lobpcg")
	eig_lobpcg = H.solve(max_states = num_states, method='lobpcg').energies
	print("Eigenvalues (lobpcg):", eig_lobpcg)

	# try to solve with some other method
	print('-'*30)
	print("Solving with foo")
	try:
		eig_foo = H.solve(max_states= num_states, method='foo')
	except NotImplementedError as e:
		print(e)

	# compare computed eigenvalues
	print('-'*30)
	print("Diff. over mean:", 2*np.abs(eig_eigsh-eig_lobpcg)/(eig_eigsh+eig_lobpcg))

