import numpy as np
from qmsolve import Hamiltonian, TwoFermions, SplitStepMethod
import qmsolve.util.constants as const


# Constants
N = 256  # Number of points to use
L = 80
X1, X2 = np.meshgrid(L*np.linspace(-0.5, 0.5 - 1.0/N, N),
                     L*np.linspace(-0.5, 0.5 - 1.0/N, N))
DX = X1[1] - X1[0]
DT = 1.0


def potential(particles):
    x1, x2 = particles.x1, particles.x2
    r = np.abs(x1 - x2)
    for i in range(N):
        r[i, i] = r[i, i - 1] if i > 0 else r[i, i+1]
    non_interacting_term = ((x1/L)**2 + (x2/L)**2)
    return const.e**2/(4.0*np.pi*const.ε0*r) + non_interacting_term


H = Hamiltonian(particles=TwoFermions(), potential=potential, N=N,
                extent=L, spatial_ndim=1)

U = SplitStepMethod(H, timestep=DT)
U.normalize_at_each_step(True)

# The wavefunction
SIGMA = 0.1
wavefunc = np.exp(-((X1/L+0.25)/SIGMA)**2/2.0
                  - ((X2/L-0.1)/SIGMA)**2/2.0
                    )*(1.0 + 0.0j) # *np.exp(20.0j*np.pi*X1/L)
wavefunc = wavefunc/np.sqrt(np.sum(wavefunc*np.conj(wavefunc)))


def main_animation():
    import matplotlib.pyplot as plt
    import matplotlib.animation as animation
    from qmsolve.util.colour_functions import complex_to_rgba
    fig = plt.figure()
    axes = fig.subplots(1, 2)
    ax = axes[0]
    ax2 = axes[1]
    max_val = np.amax(np.abs(wavefunc))
    im = ax.imshow(complex_to_rgba(X1 + 1.0j*X2),
                    extent=(X1[0, 1], X1[0, -1], X2[0, 0], X2[-1, 0]),
                    origin='lower',
                    aspect='auto')
    V = U.V
    potential_im_data = np.transpose(np.array([V, V, V, 
                                               np.amax(V)*np.ones([N, N])])
                                               /np.amax(V), (1, 2, 0))
    im2 = ax.imshow(potential_im_data,
                    origin='lower',
                    extent=(X1[0, 1], X1[0, -1], X2[0, 0], X2[-1, 0]))
    ax.set_xlabel('x1 (Hartrees)')
    ax.set_ylabel('x2 (Hartrees)')
    ax.set_title('Two 1D Particles Wavefunction')
    wavefunc_data = {'psi': (wavefunc - wavefunc.T)/np.sqrt(2.0)}
    line, = ax2.plot(X1[0], np.dot(np.abs(wavefunc)**2, np.ones([N])))
    non_interacting_term = 100.0*((X1/L)**2 + (X2/L)**2)
    ax2.plot(X1[0], np.amax(np.dot(np.abs(wavefunc)**2, np.ones([N])))*
                    non_interacting_term[N//2]/
                    np.amax(non_interacting_term[N//2]), color='gray')
    ax2.set_xlim(X1[0, 0], X1[0, -1])
    ax2.set_ylim(-0.05/10.0, 0.05*0.5)
    ax2.set_xlabel('X (Hartrees)')
    ax2.set_yticks([])

    def animation_func(*_):
        """
        Animation function
        """
        wavefunc_data['psi'] = U(wavefunc_data['psi'])
        psi = wavefunc_data['psi']
        prob = np.abs(psi)**2
        prob_1d = np.dot(prob, np.ones([N]))
        line.set_ydata(prob_1d)
        im.set_data(complex_to_rgba(psi, max_val=max_val))
        return (im, line)

    ani = animation.FuncAnimation(fig, animation_func, 
                                  blit=True, interval=1.0)
    plt.show()


main_animation()
