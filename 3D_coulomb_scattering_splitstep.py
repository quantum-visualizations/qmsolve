from qmsolve import Hamiltonian, SingleParticle, SplitStepMethod
import numpy as np
import qmsolve.util.constants as const


# Constants
N = 64  # Number of points to use
L = 4.0
S = L*np.linspace(-0.5, 0.5, N)
X, Y, Z = np.meshgrid(S, S, S)
Z_ATOMIC = 3.0 # atomic number
DT = 0.001  # timestep in seconds


# The Potential
def coulomb_potential(p):
    r = np.sqrt((p.x/L)**2 + (p.y/L)**2 + (p.z/L)**2)
    return -Z_ATOMIC*const.e**2/(4.0*np.pi*const.ε0*r)


H = Hamiltonian(particles=SingleParticle(),
                potential=coulomb_potential,
                spatial_ndim=3, N=N, extent=L)


U = SplitStepMethod(H, timestep=DT)


# The wavefunction
sigma = 0.07
# sigma = 3.0
e = (1.0 + 0.0j)*np.exp(10.0j*np.pi*X/L)
psi1 = e*np.exp(-((X/L+0.25)/sigma)**2/2.0
                -((Y/L)/sigma)**2/2.0
                -((Z/L)/sigma)**2/2.0)
psi2 = e*np.exp(-((X/L-0.25)/sigma)**2/2.0
                -((Y/L)/sigma)**2/2.0
                -((Z/L)/sigma)**2/2.0)
psi12 = psi1 + np.conj(psi2)
psi = psi12/np.sqrt(np.sum(psi12*np.conj(psi12)))
data = {'psi': 4.0*psi1/np.sqrt(np.sum(psi1*np.conj(psi1)))}


def main_animation():
    from mayavi import mlab
    mlab.figure(1, fgcolor=(1, 1, 1), bgcolor=(0, 0, 0))
    plot_data = mlab.pipeline.scalar_field(np.abs(psi))
    angle_data = (np.angle(psi)).T.ravel()
    plot_data.image_data.point_data.add_array(angle_data)
    plot_data.image_data.point_data.get_array(1).name = 'phase'
    plot_data.update()

    plot_data2 = mlab.pipeline.set_active_attribute(plot_data, point_scalars='scalar')
    contour = mlab.pipeline.contour(plot_data2)
    contour2 = mlab.pipeline.set_active_attribute(contour, point_scalars='phase')
    mlab.pipeline.surface(contour2, 
                      colormap='hsv'
                      )

    @mlab.animate(delay=10)
    def animation():
        while (1):
            for _ in range(1):
                data['psi'] = U(data['psi'])
                plot_data.mlab_source.scalars = np.abs(data['psi'])
                np.copyto(angle_data, np.angle(data['psi']).T.ravel())
                # print(np.sum(np.abs(data['psi'])**2))
            yield

    animation()
    mlab.show()


main_animation()

