import numpy as np
from qmsolve import Hamiltonian, SingleParticle, TimeSimulation, init_visualization, femtoseconds, m_e, Å

#=========================================================================================================#
#First, we define the Hamiltonian of a single particle confined in an harmonic oscillator potential. 
#=========================================================================================================#

#interaction potential
def harmonic_oscillator(particle):
    m = m_e
    T = 0.5*femtoseconds
    w = 2*np.pi/T
    k = m* w**2

    v0 = 64. * Å / femtoseconds
    print("oscillation_amplitude ", np.sqrt(m/k) *v0/Å, " amstrongs")

    return 0.5 * k * particle.x**2    +    0.5 * k * particle.y**2


#build the Hamiltonian of the system
H = Hamiltonian(particles = SingleParticle(m = m_e), 
                potential = harmonic_oscillator, 
                spatial_ndim = 2, N = 300, extent = 30 * Å)



#=========================================================================================================#
# Define the wavefunction at t = 0  (initial condition)
#=========================================================================================================#

def initial_wavefunction(particle):
    #This wavefunction correspond to a gaussian wavepacket with a mean X momentum equal to p_x0
    σ = 1.0 * Å
    v0 = 64. * Å / femtoseconds
    p_x0 = m_e * v0
    return np.exp( -1/(4* σ**2) * ((particle.x-0)**2+(particle.y-0)**2)) / np.sqrt(2*np.pi* σ**2)  *np.exp(p_x0*particle.x*1j)


#=========================================================================================================#
# Set and run the simulation
#=========================================================================================================#


total_time = 1.5 * femtoseconds
sim = TimeSimulation(hamiltonian = H, method = "split-step")
sim.run(initial_wavefunction, total_time = total_time, dt = total_time/1000., store_steps = 400)


#=========================================================================================================#
# Finally, we visualize the time dependent simulation
#=========================================================================================================#

visualization = init_visualization(sim)
visualization.animate(xlim=[-15* Å,15* Å], ylim=[-15* Å,15* Å], potential_saturation = 0.5, wavefunction_saturation = 0.2, animation_duration = 10, save_animation = False)


#for visualizing a single frame, use plot method instead of animate:
#visualization.plot(t = 5/4 * 1 * femtoseconds,xlim=[-15* Å,15* Å], ylim=[-15* Å,15* Å], potential_saturation = 0.5, wavefunction_saturation = 0.2)
