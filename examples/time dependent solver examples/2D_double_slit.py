import numpy as np
from qmsolve import Hamiltonian, SingleParticle, TimeSimulation, init_visualization, femtoseconds, m_e, Å

#=========================================================================================================#
#First, we define the Hamiltonian of a single particle
#=========================================================================================================#

#interaction potential
def double_slit(particle):
    b = 2.0* Å # slits separation
    a = 0.5* Å # slits width
    d = 0.5* Å # slits depth

    return np.where( ((particle.x < - b/2 - a) | (particle.x > b/2 + a) | ((particle.x > -b/2)  
                     & (particle.x < b/2))) & ((particle.y < d/2) & (particle.y > -d/2) ),  1e5,  0)


#build the Hamiltonian of the system
H = Hamiltonian(particles = SingleParticle(m = m_e), 
                potential = double_slit, 
                spatial_ndim = 2, N = 256, extent = 30 * Å)



#=========================================================================================================#
# Define the wavefunction at t = 0  (initial condition)
#=========================================================================================================#

def initial_wavefunction(particle):
    #This wavefunction correspond to a gaussian wavepacket with a mean Y momentum equal to p_y0
    σ = 1.0 * Å
    v0 = 80 * Å / femtoseconds
    p_y0 = m_e * v0
    return np.exp( -1/(4* σ**2) * ((particle.x-0)**2+(particle.y+8* Å)**2)) / np.sqrt(2*np.pi* σ**2)  *np.exp(p_y0*particle.y*1j)


#=========================================================================================================#
# Set and run the simulation
#=========================================================================================================#


total_time = 0.7 * femtoseconds
sim = TimeSimulation(hamiltonian = H, method = "split-step")
sim.run(initial_wavefunction, total_time = total_time, dt = total_time/8000., store_steps = 800)


#=========================================================================================================#
# Finally, we visualize the time dependent simulation
#=========================================================================================================#

visualization = init_visualization(sim)

visualization.animate(xlim=[-15* Å,15* Å], ylim=[-15* Å,15* Å], potential_saturation = 0.5, wavefunction_saturation = 0.2, animation_duration = 10, fps = 30)
