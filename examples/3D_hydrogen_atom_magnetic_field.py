import numpy as np
from qmsolve import visualization
from qmsolve import Hamiltonian, SingleParticle, Eigenstates, Å, T, eV

#==========================================================================================================================================================================
# This example computes the eigenstates of a charged particle in a magnetic field B
#
# Notes:
# The example requires to use potential_type = "matrix" in the Hamiltonian constructor, 
# which allows to use momentum operators (px and py) in the potential term.
# particle.x, particle.px , particle.y, particle.py, , particle.z, particle.pz, are treated as operators, 
# and they are discretized as matrices. Therefore, for multiplying them,
# we use the @ instead of *, because this is the operator that represents matrix multiplication.
#==========================================================================================================================================================================



#interaction potential
def coulomb_and_magnetic_interaction(particle):

    k_c = 1.0 # (e*e / (4 * np.pi * epsilon_0))  

    coulomb_interaction = - k_c * particle.r_inv
    Bx, By, Bz = [0.  ,  0.,   100.* T] 
    B_dot_L = (Bx*(-particle.py @ particle.z + particle.pz @ particle.y) +          
               By*(particle.px @ particle.z - particle.pz @ particle.x) +           
               Bz*(-particle.px @ particle.y + particle.py @ particle.x))
    
    
    𝜇 = 0.5 # e/(2*m_e)
    paramagnetic_term = -𝜇 * B_dot_L

    d = 0.125 # e**2/(8*m_e)
    diamagnetic_term = d*((-Bx*particle.y + By*particle.x)**2 + 
                          (Bx*particle.z - Bz*particle.x)**2 + 
                          (-By*particle.z + Bz*particle.y)**2)



    magnetic_interaction = diamagnetic_term  + paramagnetic_term    # shows Zeeman effect 
    return coulomb_interaction +  magnetic_interaction

H = Hamiltonian(particles=SingleParticle(),
                N=100, potential= coulomb_and_magnetic_interaction, 
                extent=30*Å, spatial_ndim=3,  potential_type = "matrix", E_min = -15*eV)




eigenstates = H.solve( max_states = 14, method ='lobpcg')


#Note: visualization with the entire complex hue isn't fully implemented yet. There is a well known bug in the red color.
v = visualization.init_visualization(eigenstates)

v.plot_type = 'contour'
#v.plot_eigenstate(12, contrast_vals = [0.02, 0.02])
v.animate(contrast_vals = [0.02, 0.02])
