import numpy as np
from qmsolve import Hamiltonian, TwoFermions, TimeSimulation, init_visualization,Å,m_e, femtoseconds


#interaction potential
def harmonic_oscillator_plus_coulomb_interaction(fermions):
	k = 0.5
	V_harmonic = 0.5*k*fermions.x1**2 + 0.5*k*fermions.x2**2 

	k = 30.83
	r = np.abs(fermions.x1 - fermions.x2)
	r = np.where(r < 0.0001, 0.0001, r)
	V_coulomb_interaction = k/ r

	return V_harmonic + V_coulomb_interaction



#build the Hamiltonian of the system
H = Hamiltonian(particles = TwoFermions(), 
				potential = harmonic_oscillator_plus_coulomb_interaction, 
				spatial_ndim = 1, N = 200, extent = 20*Å)


def initial_wavefunction(particle):
    #This wavefunction correspond to two stationary gaussian wavepackets. The wavefunction must be antisymmetric: Ψ(x1,x2) = -Ψ(x2,x1)
    σ = 0.4 * Å
    x1 = particle.x1
    x2 = particle.x2
    𝜇01 = -5.0*Å
    𝜇02 = 2.0*Å

    return (np.exp(-(x1 - 𝜇01)**2/(4*σ**2))*np.exp(-(x2 - 𝜇02)**2/(4*σ**2)) 
            - np.exp(-(x1 - 𝜇02)**2/(4*σ**2))*np.exp(-(x2 - 𝜇01)**2/(4*σ**2)))


total_time = 1. * femtoseconds
sim = TimeSimulation(hamiltonian = H, method = "split-step")
sim.run(initial_wavefunction, total_time = total_time, dt = total_time/8000., store_steps = 400)

visualization = init_visualization(sim)
visualization.plot(t = 0, xlim=[-7.5* Å,7.5* Å], potential_saturation = 0.5, wavefunction_saturation = 0.2)

visualization.animate(xlim=[-10* Å,10* Å], potential_saturation = 500, wavefunction_saturation = 0.2, animation_duration = 10, save_animation = False)

