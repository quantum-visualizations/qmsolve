import numpy as np
from qmsolve import Hamiltonian, SingleParticle, init_visualization

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
def harmonic_oscillator(particle):

    kx = 2.001 # measured in eV / (Å**2)
    ky = 2
    Bz = 2000  # measured in T

    γ = 8.794100053860817 # e/(2*m_e) * Å
    𝜇 = γ*Bz


    harmonic_interaction =  0.5 * kx * particle.x **2 +    0.5 * ky * particle.y **2
    magnetic_interaction = 𝜇* (particle.px @ particle.y - particle.py @ particle.x)
    return harmonic_interaction + magnetic_interaction



H = Hamiltonian(particles = SingleParticle(), 
                potential = harmonic_oscillator, potential_type = "matrix",
                spatial_ndim = 3, N = 400, extent = 15)


eigenstates = H.solve(max_states = 28)

print(eigenstates.energies)
visualization = init_visualization(eigenstates)
visualization.animate()
