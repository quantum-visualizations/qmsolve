"""
For scattering simulations involving the split-step method,
the time and spatial steps must be chosen carefully or else
you will encounter numerical artifacts.
"""
import numpy as np
from qmsolve import Hamiltonian, SingleParticle, SplitStepMethod


# Constants
L = 50.0
N = 256
DT = 0.08
X, Y = np.meshgrid(L*np.linspace(-0.5, 0.5 - 1.0/N, N),
                    L*np.linspace(-0.5, 0.5 - 1.0/N, N))
DX = X[1] - X[0]

#interaction potential
def double_slit_potential(particle):
    V = np.zeros([N, N])
    y0, yf = 11*N//20 - 3, 11*N//20 + 3  
    V[y0: yf, :] = 10.0
    V[y0: yf, (64-10)*N//128: (64-7)*N//128] = 0.0
    V[y0: yf, (64+7)*N//128: (64+10)*N//128] = 0.0
    return V


H = Hamiltonian(particles = SingleParticle(), 
				potential = double_slit_potential, 
				spatial_ndim = 2, N = N, extent = L)


U = SplitStepMethod(H, timestep=DT)
# To advance an initial wavefunction psi_0:
# psi_1 = U(psi_0)

# The wavefunction
SIGMA = 0.068
wavefunc = np.exp(-((X/L)/SIGMA)**2/2.0
                    - ((Y/L-0.25)/SIGMA)**2/2.0)*np.exp(-50.0j*np.pi*Y/L)
wavefunc = wavefunc/np.sqrt(np.sum(wavefunc*np.conj(wavefunc)))


def animation():

    import matplotlib.pyplot as plt
    import matplotlib.animation as animation
    from qmsolve.util.colour_functions import complex_to_rgb

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    # print(np.amax(np.angle(psi)))
    max_val = np.amax(np.abs(wavefunc))
    im = ax.imshow(complex_to_rgb(X + 1.0j*Y),
                    extent=(X[0, 1], X[0, -1], Y[0, 0], Y[-1, 0]))
    V = U.V.T
    potential_im_data = np.transpose(np.array([V/np.amax(V), 
                                            V/np.amax(V), V/np.amax(V), 
                                            V/np.amax(V)]), (2, 1, 0))
    im2 = ax.imshow(potential_im_data,
                    extent=(X[0, 1], X[0, -1], Y[0, 0], Y[-1, 0]))
    ax.set_xlabel('x (Hartrees)')
    ax.set_ylabel('y (Hartrees)')
    ax.set_title('Wavefunction')
    data = {'psi': wavefunc}

    def animation_func(*_):
        """
        Animation function
        """
        data['psi'] = U(data['psi'])
        im.set_data(complex_to_rgb(data['psi']))
        return (im, im2)

    ani = animation.FuncAnimation(fig, animation_func, blit=True, interval=1.0)
    plt.show()


animation()
