import numpy as np
from qmsolve import Hamiltonian, SingleParticle,TimeSimulation, init_visualization, m_e, e, Å, T , eV, femtoseconds

#==========================================================================================================================================================================
# This example simulated the time evolution of a charged particle in a magnetic field Bz, that points in the z direction.
#
# Notes:
# The example requires to use potential_type = "matrix" in the Hamiltonian constructor, 
# which allows to use momentum operators (px and py) in the potential term
# particle.x, particle.px , particle.y, particle.py, are treated as operators, and they are discretized as matrices. Therefore, for multiplying them,
# we use the @ instead of *, because this is the operator that represents matrix multiplication.
#==========================================================================================================================================================================


#interaction potential
def constant_magnetic_field(particle):

    Bz = 100000 * T 

    B_dot_L =  Bz*(-particle.px @ particle.y + particle.py @ particle.x)
    𝜇 = 0.5 # e/(2*m_e)
    paramagnetic_term = -𝜇 * B_dot_L

    d = 0.125 # e**2/(8*m_e)
    diamagnetic_term = d* Bz**2 *(particle.x**2 + particle.y**2)


    magnetic_interaction = diamagnetic_term  + paramagnetic_term

    v0 = 80 * Å / femtoseconds #mean initial velocity of the wavepacket
    R = (m_e*v0 / (e*Bz))/Å #cyclotron radius
    print( (u"classical cyclotron radius = {} angstroms".format("%.2f"  % (R)) )  )

    P = (2*np.pi*m_e/(e*Bz))/ femtoseconds #cyclotron period
    print( (u"classical cyclotron period = {} femtoseconds".format("%.2f"  % (P)) )  )


    return magnetic_interaction 



H = Hamiltonian(particles = SingleParticle(), 
                potential = constant_magnetic_field, potential_type = "matrix",
                spatial_ndim = 2, N = 256, extent = 35 * Å)

#=========================================================================================================#
# Define the wavefunction at t = 0  (initial condition)
#=========================================================================================================#

def initial_wavefunction(particle):
    #This wavefunction correspond to a gaussian wavepacket with a mean X momentum equal to p_x0
    σ = 1.0 * Å
    v0 = 80 * Å / femtoseconds
    p_x0 = m_e * v0
    return np.exp( -1/(4* σ**2) * ((particle.x-0)**2+(particle.y-5.0*Å)**2)) / np.sqrt(2*np.pi* σ**2)  *np.exp(p_x0*particle.x*1j)


#=========================================================================================================#
# Set and run the simulation
#=========================================================================================================#


total_time = 1.0 * femtoseconds 
sim = TimeSimulation(hamiltonian = H, method = "crank-nicolson")
sim.run(initial_wavefunction, total_time = total_time, dt = total_time/(4000), store_steps = 200)


#=========================================================================================================#
# Finally, we visualize the time dependent simulation
#=========================================================================================================#

visualization = init_visualization(sim)


visualization.animate(xlim=[-15* Å,15* Å], ylim=[-15* Å,15* Å], potential_saturation = 0.5, wavefunction_saturation = 0.2, animation_duration = 10, fps = 30, save_animation = False)


#for visualizing a single frame, use plot method instead of animate:
#visualization.plot(t = 0.36 * femtoseconds,xlim=[-15* Å,15* Å], ylim=[-15* Å,15* Å], potential_saturation = 0.5, wavefunction_saturation = 0.2)
