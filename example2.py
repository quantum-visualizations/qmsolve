import numpy as np
from qmsolve import Halmitonian, Two_fermions


#interaction potential
def testing(fermions):

	V = fermions.x1 + fermions.x2 
	return V



H = Halmitonian(particles = Two_fermions(), 
				potential = testing, 
				spatial_ndim = 1, N = 60, extent = 10)


energies, eigenstates = H.solve(max_states = 30)

#Not tested yet. Work in progress

