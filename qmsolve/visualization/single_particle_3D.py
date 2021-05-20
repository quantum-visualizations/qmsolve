import numpy as np
from mayavi import mlab
from .visualization import Visualization


class VisualizationSingleParticle3D(Visualization):

    def __init__(self,eigenstates):
        self.eigenstates = eigenstates
        self.plot_type = 'volume'

    def slider_plot(self):
        raise NotImplementedError

    def plot_eigenstate(self, k, contrast_vals= [0.1, 0.25]):
        eigenstates = self.eigenstates.array
        mlab.figure(1, bgcolor=(0, 0, 0), size=(700, 700))
        psi = eigenstates[k]
        min_ = psi.min()
        max_ = psi.max()


        if self.plot_type == 'volume':
            
            max_= np.amax(eigenstates)
            psi = (psi)/(max_)

            L = self.eigenstates.extent/2
            N = self.eigenstates.N

            vol = mlab.pipeline.volume(mlab.pipeline.scalar_field(psi))

            # Change the color transfer function
            from tvtk.util import ctf
            c = ctf.save_ctfs(vol._volume_property)
            c['rgb'] = [[-0.45, 0.3, 0.3, 1.0],
                        [-0.4, 0.1, 0.1, 1.0],
                        [-0.3, 0.0, 0.0, 1.0],
                        [-0.2, 0.0, 0.0, 1.0],
                        [-0.001, 0.0, 0.0, 1.0],
                        [0.0, 0.0, 0.0, 0.0],
                        [0.001, 1.0, 0.0, 0.],
                        [0.2, 1.0, 0.0, 0.0],
                        [0.3, 1.0, 0.0, 0.0],
                        [0.4, 1.0, 0.1, 0.1],
                        [0.45, 1.0, 0.3, 0.3]]

            c['alpha'] = [[-0.5, 1.0],
                          [-contrast_vals[1], 1.0],
                          [-contrast_vals[0], 0.0],
                          [0, 0.0],
                          [contrast_vals[0], 0.0],
                          [contrast_vals[1], 1.0],
                         [0.5, 1.0]]
            ctf.load_ctfs(c, vol._volume_property)
            # Update the shadow LUT of the volume module.
            vol.update_ctf = True

            mlab.outline()
            mlab.axes(xlabel='x [Å]', ylabel='y [Å]', zlabel='z [Å]',nb_labels=6 , ranges = (-L,L,-L,L,-L,L) )
            #azimuth angle
            φ = 30
            mlab.view(azimuth= φ,  distance=N*3.5)
            mlab.show()

        else:

            # Adding colour scale to a contour plot:
            # https://github.com/enthought/mayavi/
            #       blob/master/examples/mayavi/mlab/atomic_orbital.py
            field = mlab.pipeline.scalar_field(np.abs(psi), 
                                            vmin=min_ + 0.008 * (max_ - min_),
                                            vmax=min_ + 0.5 * (max_ - min_))
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

    def animate(self,  contrast_vals= [0.1, 0.25]):
        eigenstates = self.eigenstates.array
        energies = self.eigenstates.energies
        mlab.figure(1, bgcolor=(0, 0, 0), size=(700, 700))

        
        if self.plot_type == 'volume':
            psi = eigenstates[0]
            max_ = np.amax(eigenstates)

            psi = (psi)/(max_)

            L = self.eigenstates.extent/2
            N = self.eigenstates.N
            field = mlab.pipeline.scalar_field(psi)
            vol = mlab.pipeline.volume(field)

            # Change the color transfer function
            from tvtk.util import ctf
            c = ctf.save_ctfs(vol._volume_property)
            c['rgb'] = [[-0.45, 0.3, 0.3, 1.0],
                        [-0.4, 0.1, 0.1, 1.0],
                        [-0.3, 0.0, 0.0, 1.0],
                        [-0.2, 0.0, 0.0, 1.0],
                        [-0.001, 0.0, 0.0, 1.0],
                        [0.0, 0.0, 0.0, 0.0],
                        [0.001, 1.0, 0.0, 0.],
                        [0.2, 1.0, 0.0, 0.0],
                        [0.3, 1.0, 0.0, 0.0],
                        [0.4, 1.0, 0.1, 0.1],
                        [0.45, 1.0, 0.3, 0.3]]

            c['alpha'] = [[-0.5, 1.0],
                          [-contrast_vals[1], 1.0],
                          [-contrast_vals[0], 0.0],
                          [0, 0.0],
                          [contrast_vals[0], 0.0],
                          [contrast_vals[1], 1.0],
                         [0.5, 1.0]]
            ctf.load_ctfs(c, vol._volume_property)
            # Update the shadow LUT of the volume module.
            vol.update_ctf = True

            mlab.outline()
            mlab.axes(xlabel='x [Å]', ylabel='y [Å]', zlabel='z [Å]',nb_labels=6 , ranges = (-L,L,-L,L,-L,L) )

            #azimuth angle
            φ = 30
            mlab.view(azimuth= φ,  distance=N*3.5)


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

                    psi = (psi)/(max_)
                    field.mlab_source.scalars = psi
                    # Change the color transfer function
                    from tvtk.util import ctf
                    c = ctf.save_ctfs(vol._volume_property)
                    c['rgb'] = [[-0.45, 0.3, 0.3, 1.0],
                                [-0.4, 0.1, 0.1, 1.0],
                                [-0.3, 0.0, 0.0, 1.0],
                                [-0.2, 0.0, 0.0, 1.0],
                                [-0.001, 0.0, 0.0, 1.0],
                                [0.0, 0.0, 0.0, 0.0],
                                [0.001, 1.0, 0.0, 0.],
                                [0.2, 1.0, 0.0, 0.0],
                                [0.3, 1.0, 0.0, 0.0],
                                [0.4, 1.0, 0.1, 0.1],
                                [0.45, 1.0, 0.3, 0.3]]

                    c['alpha'] = [[-0.5, 1.0],
                                  [-contrast_vals[1], 1.0],
                                  [-contrast_vals[0], 0.0],
                                  [0, 0.0],
                                  [contrast_vals[0], 0.0],
                                  [contrast_vals[1], 1.0],
                                 [0.5, 1.0]]
                    ctf.load_ctfs(c, vol._volume_property)
                    # Update the shadow LUT of the volume module.
                    vol.update_ctf = True

                    φ = 30 + data['t'] * 360 / 10 
                    mlab.view(azimuth= φ, distance=N*3.5)

                    yield

            animation()
            mlab.show()


        elif self.plot_type == 'volume_old_colormap':
            L = self.eigenstates.extent/2
            N = self.eigenstates.N

            xx, yy, zz = np.mgrid[-L:L:N*1j, -L:L:N*1j, -L:L:N*1j]
            psi = np.abs(eigenstates[0])
            min_ = np.amin(np.abs(eigenstates[:]))
            max_ = np.amax(np.abs(eigenstates[:]))

            vmin = min_ + 0.008 * (max_ - min_)
            vmax = min_ + 0.4 * (max_ - min_)
            psi = np.where(psi > vmax, vmax,psi)
            psi = np.where(psi < vmin, vmin,psi)

            field = mlab.pipeline.scalar_field(xx, yy, zz, psi)
            vol = mlab.pipeline.volume(field)
            φ = 30
            mlab.view(azimuth= φ, distance=L*7.2)
            mlab.outline()

            mlab.axes(xlabel='x [Å]', ylabel='y [Å]', zlabel='z [Å]',nb_labels=6 , ranges = (-L,L,-L,L,-L,L) )

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
                    vmin = min_ + 0.008 * (max_ - min_)
                    vmax = min_ + 0.4 * (max_ - min_)
                    psi = np.where(psi > vmax, vmax,psi)
                    psi = np.where(psi < vmin, vmin,psi)

                    field.mlab_source.scalars = np.abs(psi)
                    φ = 30 + data['t'] * 360 / 10 
                    mlab.view(azimuth= φ, distance=L*7.2)

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

    def superpositions(self, states, **kw):
        params = {'contrast_vals': [0.1, 0.25], 'dt': 0.1,
                  'display_type': 'real'}
        for k in kw.keys():
            if k in params:
                params[k] = kw[k]
            else:
                raise KeyError
        eigenstates = self.eigenstates.array[0 : len(states)]
        coeffs = states
        energies = self.eigenstates.energies
        mlab.figure(1, bgcolor=(0, 0, 0), size=(700, 700))
        psi = sum([eigenstates[i]*coeffs[i]
                   for i in range(len(eigenstates))])
        max_ = np.amax(np.abs(psi))
        if params['display_type'] == 'real':
            psi = np.real(psi)/(max_)
        elif params['display_type'] == 'imag':
            psi = np.imag(psi)/(max_)
        else:
            psi = np.abs(psi)/(max_)

        L = self.eigenstates.extent/2
        N = self.eigenstates.N
        field = mlab.pipeline.scalar_field(psi)
        vol = mlab.pipeline.volume(field)

        # Change the color transfer function
        from tvtk.util import ctf
        c = ctf.save_ctfs(vol._volume_property)
        c['rgb'] = [[-0.45, 0.3, 0.3, 1.0],
                    [-0.4, 0.1, 0.1, 1.0],
                    [-0.3, 0.0, 0.0, 1.0],
                    [-0.2, 0.0, 0.0, 1.0],
                    [-0.001, 0.0, 0.0, 1.0],
                    [0.0, 0.0, 0.0, 0.0],
                    [0.001, 1.0, 0.0, 0.],
                    [0.2, 1.0, 0.0, 0.0],
                    [0.3, 1.0, 0.0, 0.0],
                    [0.4, 1.0, 0.1, 0.1],
                    [0.45, 1.0, 0.3, 0.3]]

        contrast_vals = params['contrast_vals']
        c['alpha'] = [[-0.5, 1.0],
                        [-contrast_vals[1], 1.0],
                        [-contrast_vals[0], 0.0],
                        [0, 0.0],
                        [contrast_vals[0], 0.0],
                        [contrast_vals[1], 1.0],
                        [0.5, 1.0]]
        ctf.load_ctfs(c, vol._volume_property)
        # Update the shadow LUT of the volume module.
        vol.update_ctf = True

        mlab.outline()
        mlab.axes(xlabel='x [Å]', ylabel='y [Å]', zlabel='z [Å]',nb_labels=6 , 
                  ranges = (-L,L,-L,L,-L,L) )

        #azimuth angle
        φ = 30
        mlab.view(azimuth= φ,  distance=N*3.5)


        data = {'t': 0.0}
        @mlab.animate(delay=10)
        def animation():
            while (1):
                data['t'] += params['dt']
                t = data['t']
                psi = sum([eigenstates[i]*np.exp(-1.0j*energies[i]*t)*coeffs[i]
                           for i in range(len(eigenstates))])
                if params['display_type'] == 'real':
                    psi = np.real(psi)/(max_)
                elif params['display_type'] == 'imag':
                    psi = np.imag(psi)/(max_)
                else:
                    psi = np.abs(psi)/(max_)
                field.mlab_source.scalars = psi
                # Change the color transfer function
                from tvtk.util import ctf
                c = ctf.save_ctfs(vol._volume_property)
                c['rgb'] = [[-0.45, 0.3, 0.3, 1.0],
                            [-0.4, 0.1, 0.1, 1.0],
                            [-0.3, 0.0, 0.0, 1.0],
                            [-0.2, 0.0, 0.0, 1.0],
                            [-0.001, 0.0, 0.0, 1.0],
                            [0.0, 0.0, 0.0, 0.0],
                            [0.001, 1.0, 0.0, 0.],
                            [0.2, 1.0, 0.0, 0.0],
                            [0.3, 1.0, 0.0, 0.0],
                            [0.4, 1.0, 0.1, 0.1],
                            [0.45, 1.0, 0.3, 0.3]]

                c['alpha'] = [[-0.5, 1.0],
                                [-contrast_vals[1], 1.0],
                                [-contrast_vals[0], 0.0],
                                [0, 0.0],
                                [contrast_vals[0], 0.0],
                                [contrast_vals[1], 1.0],
                                [0.5, 1.0]]
                ctf.load_ctfs(c, vol._volume_property)
                # Update the shadow LUT of the volume module.
                vol.update_ctf = True

                φ = 30 + data['t'] * 360 / 10 
                mlab.view(azimuth= φ, distance=N*3.5)

                yield

        animation()
        mlab.show()

