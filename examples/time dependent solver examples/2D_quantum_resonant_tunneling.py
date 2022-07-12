import numpy as np
from qmsolve import Hamiltonian, SingleParticle, TimeSimulation, init_visualization, femtoseconds, m_e, Å, eV

#=========================================================================================================#
# First, we define the Hamiltonian of a single particle. 
# The interaction potential consists of a double potential barrier
#=========================================================================================================#


def potential(particle):
    
    #interaction potential
    a = 3*Å
    b = 5*Å
    V0 = 1.7*eV
    
    #first barrier
    V =  np.where((particle.x>0) & (particle.x<a) , V0, 0.)
    
    #second barrier
    V =  np.where((particle.x>a+b) & (particle.x<a+b+a)  , V + V0 , V)
    return V


#build the Hamiltonian of the system
H = Hamiltonian(particles = SingleParticle(m = m_e), 
                potential = potential, 
                spatial_ndim = 2, N = 256, extent = 700*Å)



#=========================================================================================================#
# Define the wavefunction at t = 0  (initial condition)
#=========================================================================================================#

E1 = 0.6*eV
E2 = 1.0*eV

σ = 25.0*Å
p1_x0 = np.sqrt(2*E1 / m_e)
p2_x0 = np.sqrt(2*E2 / m_e)

def initial_wavefunction(particle):
    #This wavefunction correspond to a gaussian wavepacket with a mean X momentum equal to p_x0
    return (np.exp( -1/(4* σ**2) * ((particle.x+100*Å)**2+ 1*(particle.y-250)**2)) / np.sqrt(2*np.pi* σ**2)  *np.exp(p1_x0*particle.x*1j) +
           np.exp( -1/(4* σ**2) * ((particle.x+100*Å)**2+ 1*(particle.y+250)**2)) / np.sqrt(2*np.pi* σ**2)  *np.exp(p2_x0*particle.x*1j))


#=========================================================================================================#
# Set and run the simulation
#=========================================================================================================#


total_time = 4000
sim = TimeSimulation(hamiltonian = H, method = "split-step")
sim.run(initial_wavefunction, total_time = total_time, dt = total_time/8000., store_steps = 800)


#=========================================================================================================#
# Finally, we visualize the time dependent simulation
#=========================================================================================================#

visualization = init_visualization(sim)

visualization.animate(xlim=[-300*Å,300*Å], ylim=[-300*Å,300*Å], potential_saturation = 0.3, wavefunction_saturation = 0.3, animation_duration = 10, fps = 30)
