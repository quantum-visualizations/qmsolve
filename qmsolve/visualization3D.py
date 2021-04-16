import numpy as np
from mayavi import mlab


def visualize3D(energies, eigenstates, k, plot_type='volume'):
    mlab.figure(1, bgcolor=(0, 0, 0), size=(350, 350))

    psi2 = np.abs(eigenstates[k])
    min = psi2.min()
    max = psi2.max()
    if plot_type == 'volume':

        vol = mlab.pipeline.volume(mlab.pipeline.scalar_field(psi2), 
                                vmin=min + 0.008 * (max - min),        
                                vmax=min + 0.5 * (max - min))

    else:

        # Adding colour scale to a contour plot:
        # https://github.com/enthought/mayavi/
        #       blob/master/examples/mayavi/mlab/atomic_orbital.py
        field = mlab.pipeline.scalar_field(psi2, 
                                           vmin=min + 0.008 * (max - min),
                                           vmax=min + 0.5 * (max - min))
        field.image_data.point_data.add_array(eigenstates[k].T.ravel())
        field.image_data.point_data.get_array(1).name = 'phase'
        field.update()
        field2 = mlab.pipeline.set_active_attribute(field, 
                                                    point_scalars='scalar')
        contour = mlab.pipeline.contour(field2)
        contour2 = mlab.pipeline.set_active_attribute(contour, 
                                                      point_scalars='phase')
        mlab.pipeline.surface(contour2)

    mlab.outline()
    mlab.axes(nb_labels=4)


    φ = 30
    mlab.view(azimuth= φ)

    mlab.show()


def animate3D(energies, eigenstates):
    mlab.figure(1, bgcolor=(0, 0, 0), size=(350, 350))
    field = mlab.pipeline.scalar_field(np.abs(eigenstates[0]))
    colour_data = np.angle(eigenstates[0].T.ravel())
    field.image_data.point_data.add_array(colour_data)
    field.image_data.point_data.get_array(1).name = 'phase'
    field.update()
    field2 = mlab.pipeline.set_active_attribute(field, 
                                                point_scalars='scalar')
    contour = mlab.pipeline.contour(field2)
    contour2 = mlab.pipeline.set_active_attribute(contour, 
                                                  point_scalars='phase')
    mlab.pipeline.surface(contour2)
    data = {'t': 0.0}
    @mlab.animate(delay=10)
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
            np.copyto(colour_data, np.angle(psi.T.ravel()))
            field.mlab_source.scalars = np.abs(psi)
            yield
    animation()
    mlab.show()
