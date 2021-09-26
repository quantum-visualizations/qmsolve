import numpy as np
from qmsolve import Hamiltonian, SingleParticle, TimeSimulation, init_visualization, femtoseconds, m_e, Å


def lattice_point(x,y, x0,y0):
    σ = 0.5
    V0 = 40
    return V0*np.exp( -1/(4* σ**2) * ((x-x0)**2+(y-y0)**2))

Nx_point = 12
Ny_point = 12

x_point, y_point = np.meshgrid(np.linspace(-50/2, 50/2, Nx_point),
                               np.linspace(-50/2, 50/2, Ny_point))

#interaction potential
def lattice(particle):
    V = 0
    for i in range(0,Nx_point):
        for j in range(Ny_point//2,Ny_point):
            V += lattice_point(particle.x,particle.y, x_point[i,j],y_point[i,j])
    return V


#build the Hamiltonian of the system
H = Hamiltonian(particles = SingleParticle(m = m_e), 
                potential = lattice, 
                spatial_ndim = 2, N = 300, extent = 50 * Å)


# Define the wavefunction at t = 0  (initial condition)
def initial_wavefunction(particle):
    #This wavefunction correspond to a gaussian wavepacket with a mean X momentum equal to p_x0
    σ = 1.0 * Å
    v0 = 80 * Å / femtoseconds
    p_x0 = m_e * v0
    return np.exp( -1/(4* σ**2) * ((particle.x+15)**2+(particle.y-0)**2)) / np.sqrt(2*np.pi* σ**2)  *np.exp(p_x0*particle.x*1j)



# Set and run the simulation
total_time = 0.4 * femtoseconds
sim = TimeSimulation(hamiltonian = H, method = "split-step")
sim.run(initial_wavefunction, total_time = total_time, dt = total_time/2000., store_steps = 500)


# Finally, we visualize the time dependent simulation
visualization = init_visualization(sim)
visualization.animate(xlim=[-25* Å,25* Å], ylim=[-25* Å,25* Å], potential_saturation = 1.0, wavefunction_saturation = 0.1, animation_duration = 10, fps = 20)
