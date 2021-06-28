import numpy as np
from qmsolve import Hamiltonian, SingleParticle, SplitStepMethod


L = 1.0
N = 256
DT = 0.0005
X, Y = np.meshgrid(L*np.linspace(-0.5, 0.5 - 1.0/N, N),
                    L*np.linspace(-0.5, 0.5 - 1.0/N, N))
DX = X[1] - X[0]

#interaction potential
def harmonic_oscillator(particle):

	kx =  31250.0
	ky =  31250.0
	return 0.5 * kx * particle.x**2    +    0.5 * ky * particle.y**2

H = Hamiltonian(particles = SingleParticle(), 
				potential = harmonic_oscillator, 
				spatial_ndim = 2, N = N, extent = L)


U = SplitStepMethod(H, timestep=0.0005)
# To advance an initial wavefunction psi_0:
# psi_1 = U(psi_0)
# psi must have the same dimensions as the potential.

# The wavefunction
SIGMA = 0.068
wavefunc = np.exp(-((X/L+0.25)/SIGMA)**2/2.0
                    - ((Y/L-0.25)/SIGMA)**2/2.0)*(1.0 + 0.0j)
wavefunc = wavefunc/np.sqrt(np.sum(wavefunc*np.conj(wavefunc)))

def animation():

    import matplotlib.pyplot as plt
    import matplotlib.animation as animation
    from qmsolve.util.colour_functions import complex_to_rgba

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    # print(np.amax(np.angle(psi)))
    max_val = np.amax(np.abs(wavefunc))
    im = ax.imshow(complex_to_rgba(X + 1.0j*Y),
                    extent=(X[0, 1], X[0, -1], Y[0, 0], Y[-1, 0]))
    V = U.V
    potential_im_data = np.transpose(np.array([V, V, V, np.amax(V)*np.ones([N, N])])
                                        /np.amax(V), (2, 1, 0))
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
        im.set_data(complex_to_rgba(data['psi'], max_val=max_val))
        return (im,)

    ani = animation.FuncAnimation(fig, animation_func, blit=True, interval=1.0)
    plt.show()

animation()

