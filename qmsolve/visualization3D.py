import numpy as np
from mayavi import mlab


def visualize3D(energies, eigenstates, k):
    mlab.figure(1, bgcolor=(0, 0, 0), size=(350, 350))

    psi2 = np.abs(eigenstates[k])
    min = psi2.min()
    max = psi2.max()

    vol = mlab.pipeline.volume(mlab.pipeline.scalar_field(psi2), 
                               vmin=min + 0.008 * (max - min),        
                               vmax=min + 0.5 * (max - min))

    
    mlab.outline()
    mlab.axes(nb_labels=4)


    φ = 30
    mlab.view(azimuth= φ)

    mlab.show()



def animate3D(energies, eigenstates):
    # mlab.figure(1, fgcolor=(1, 1, 1), bgcolor=(0, 0, 0))
    plot_data = mlab.contour3d(np.abs(eigenstates[0]))
    data = {'t': 0.0}
    @mlab.animate
    def animation():
        while (1):
            data['t'] += 0.05
            k1 = int(data['t']) % len(energies)
            k2 = (int(data['t']) + 1) % len(energies)
            if data['t'] % 1.0 > 0.5:
                t = (data['t'] - int(data['t']) - 0.5)
                psi = (np.cos(np.pi*t)*eigenstates[k1]
                       + np.sin(np.pi*t)*eigenstates[k2])
            else:
                psi = eigenstates[k1]
            plot_data.mlab_source.scalars = np.abs(psi)
            yield
    animation()
    mlab.show()
