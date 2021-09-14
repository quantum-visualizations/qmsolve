import numpy as np
from qmsolve import Hamiltonian, SingleParticle, TimeSimulation, init_visualization, femtoseconds, m_e, Å


#interaction potential
def harmonic_oscillator(particle):
    m = m_e
    T = 5*femtoseconds
    w = 2*np.pi/T
    k = m* w**2
    return 0.5 * k * particle.x**2    +    0.5 * k * particle.y**2


#build the Hamiltonian of the system
H = Hamiltonian(particles = SingleParticle(m = m_e), 
                potential = harmonic_oscillator, 
                spatial_ndim = 2, N = 300, extent = 30 * Å)


#wavefunction at t = 0. 
def initial_wavefunction(particle):
    #This wavefunction correspond to a gaussian wavepacket with a mean X momentum equal to p_x0
    σ = 1.0 * Å
    v0 = 8 * Å / femtoseconds
    p_x0 = m_e * v0
    return np.exp( -1/(4* σ**2) * ((particle.x-0)**2+(particle.y-0)**2)) / np.sqrt(2*np.pi* σ**2)  *np.exp(p_x0*particle.x*1j)


total_time = 15 * femtoseconds

#set the time dependent simulation
sim = TimeSimulation(hamiltonian = H, method = "split-step-cupy")
sim.run(initial_wavefunction, total_time = total_time, dt = total_time/1000., store_steps = 400)

#visualize the time dependent simulation
visualization = init_visualization(sim)
visualization.animate(xlim=[-15* Å,15* Å], ylim=[-15* Å,15* Å], potential_saturation = 0.5, wavefunction_saturation = 0.2, animation_duration = 10, save_animation = False)


#for visualizing a single frame, use plot method instead of animate:
#visualization.plot(t = 5/4 * 1 * femtoseconds,xlim=[-15* Å,15* Å], ylim=[-15* Å,15* Å], potential_saturation = 0.5, wavefunction_saturation = 0.2)
