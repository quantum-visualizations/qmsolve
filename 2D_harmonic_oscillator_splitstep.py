import numpy as np
from qmsolve import Hamiltonian, SingleParticle, SplitStepMethod



#interaction potential
def harmonic_oscillator(particle):

	kx = 2 # measured in eV / (Å**2)
	ky = 2
	return 0.5 * kx * particle.x**2    +    0.5 * ky * particle.y**2


L = 2.0
N = 200

H = Hamiltonian(particles = SingleParticle(), 
				potential = harmonic_oscillator,
				spatial_ndim = 2, N = N, extent = L)


U = SplitStepMethod(H, 1e-10)
V = U.V

import matplotlib.pyplot as plt
import matplotlib.animation as animation

X, Y = np.meshgrid(L*np.linspace(-0.5, 0.5 - 1.0/N, N),
                    L*np.linspace(-0.5, 0.5 - 1.0/N, N))
DX = X[1] - X[0]

# The wavefunction
SIGMA = 0.056568
wavefunc = np.exp(-((X/L+0.25)/SIGMA)**2/2.0
                    - ((Y/L-0.25)/SIGMA)**2/2.0)*(1.0 + 0.0j)
wavefunc = wavefunc/np.sqrt(np.sum(wavefunc*np.conj(wavefunc)))

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
# print(np.amax(np.angle(psi)))
max_val = np.amax(np.abs(wavefunc))
im = ax.imshow(np.angle(X + 1.0j*Y),
                alpha=np.abs(wavefunc)/max_val,
                extent=(X[0, 1], X[0, -1], Y[0, 0], Y[-1, 0]),
                interpolation='none',
                cmap='hsv')
potential_im_data = np.transpose(np.array([V, V, V, np.amax(V)*np.ones([N, N])])
                                    /np.amax(V), (2, 1, 0))
im2 = ax.imshow(potential_im_data,
                extent=(X[0, 1], X[0, -1], Y[0, 0], Y[-1, 0]),
                interpolation='bilinear')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_title('Wavefunction')
data = {'psi': wavefunc, 'steps': 0}

def animation_func(*_):
    """
    Animation function
    """
    data['psi'] = U(data['psi'])
    im.set_data(np.angle(data['psi']))
    im.set_alpha(np.abs(data['psi'])/max_val)
    data['steps'] += 1
    return (im2, im)


ani = animation.FuncAnimation(fig, animation_func, blit=True, interval=1.0)
plt.show()



