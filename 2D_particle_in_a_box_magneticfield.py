import numpy as np
from qmsolve import Hamiltonian, SingleParticle, init_visualization, Å, T , eV

#==========================================================================================================================================================================
# This example computes the eigenstates of a charged particle in a magnetic field Bz, that points in the z direction.
#
# Notes:
# The example requires to use potential_type = "matrix" in the Hamiltonian constructor, 
# which allows to use momentum operators (px and py) in the potential term
# particle.x, particle.px , particle.y, particle.py, are treated as operators, and they are discretized as matrices. Therefore, for multiplying them,
# we use the @ instead of *, because this is the operator that represents matrix multiplication.
#==========================================================================================================================================================================


#interaction potential
def magnetic_interaction(particle):

    Bz = 8000 * T 

    γ = 0.5 # e/(2*m_e)
    𝜇 = γ*Bz

    magnetic_interaction = 𝜇* (particle.px @ particle.y - particle.py @ particle.x)
    return particle.x*0.0000001 + magnetic_interaction



H = Hamiltonian(particles = SingleParticle(), 
                potential = magnetic_interaction, potential_type = "matrix",
                spatial_ndim = 2, N = 400, extent = 15 * Å, E_min=0.0)


eigenstates = H.solve(max_states = 28)

print(eigenstates.energies)
visualization = init_visualization(eigenstates)
visualization.animate()
