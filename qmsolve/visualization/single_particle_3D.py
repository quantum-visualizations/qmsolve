import numpy as np
from mayavi import mlab
from .visualization import Visualization


class VisualizationSingleParticle3D(Visualization):

    def __init__(self,eigenstates):
        self.eigenstates = eigenstates
        self.plot_type = 'volume'

    def slider_plot(self):
        raise NotImplementedError

    def plot_eigenstate(self, k):
        eigenstates = self.eigenstates.array
        mlab.figure(1, bgcolor=(0, 0, 0), size=(350, 350))
        psi = np.abs(eigenstates[k])
        min = psi.min()
        max = psi.max()



        if self.plot_type == 'volume':

            L = self.eigenstates.extent/2
            N = self.eigenstates.N

            xx, yy, zz = np.mgrid[-L:L:N*1j, -L:L:N*1j, -L:L:N*1j]
            vol = mlab.pipeline.volume(mlab.pipeline.scalar_field(xx,yy,zz,psi), 
                                    vmin=min + 0.008 * (max - min),        
                                    vmax=min + 0.5 * (max - min))

            mlab.outline()

            mlab.axes(xlabel='x [Å]', ylabel='y [Å]', zlabel='z [Å]',nb_labels=4)

        else:

            # Adding colour scale to a contour plot:
            # https://github.com/enthought/mayavi/
            #       blob/master/examples/mayavi/mlab/atomic_orbital.py
            field = mlab.pipeline.scalar_field(psi2, 
                                            vmin=min + 0.008 * (max - min),
                                            vmax=min + 0.5 * (max - min))
            field.image_data.point_data.add_array(eigenstates[k].ravel())
            field.image_data.point_data.get_array(1).name = 'phase'
            field.update()
            field2 = mlab.pipeline.set_active_attribute(field, 
                                                        point_scalars='scalar')
            contour = mlab.pipeline.contour(field2)
            contour2 = mlab.pipeline.set_active_attribute(contour, 
                                                        point_scalars='phase')
            mlab.pipeline.surface(contour2)

        mlab.show()

    def animate(self):
        eigenstates = self.eigenstates.array
        energies = self.eigenstates.energies
        mlab.figure(1, bgcolor=(0, 0, 0), size=(350, 350))


        if self.plot_type == 'volume':
            L = self.eigenstates.extent/2
            N = self.eigenstates.N

            xx, yy, zz = np.mgrid[-L:L:N*1j, -L:L:N*1j, -L:L:N*1j]
            psi = np.abs(eigenstates[0])
            min = np.amin(np.abs(eigenstates[:]))
            max = np.amax(np.abs(eigenstates[:]))

            vmin = min + 0.008 * (max - min)
            vmax = min + 0.4 * (max - min)
            psi = np.where(psi > vmax, vmax,psi)
            psi = np.where(psi < vmin, vmin,psi)

            field = mlab.pipeline.scalar_field(xx, yy, zz, psi)
            vol = mlab.pipeline.volume(field)
            φ = 30
            mlab.view(azimuth= φ, distance=L*6.7)
            mlab.outline()

            mlab.axes(xlabel='x [Å]', ylabel='y [Å]', zlabel='z [Å]',nb_labels=4)

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

                    psi = np.abs(psi)
                    vmin = min + 0.008 * (max - min)
                    vmax = min + 0.4 * (max - min)
                    psi = np.where(psi > vmax, vmax,psi)
                    psi = np.where(psi < vmin, vmin,psi)

                    field.mlab_source.scalars = np.abs(psi)
                    φ = 30 + data['t'] * 360 / 10 
                    mlab.view(azimuth= φ, distance=L*6.7)

                    yield

            animation()
            mlab.show()

        else:

            field = mlab.pipeline.scalar_field(np.abs(eigenstates[0]))
            colour_data = np.angle(eigenstates[0].ravel())
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
