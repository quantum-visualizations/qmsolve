import numpy as np
from mayavi import mlab
from .visualization import Visualization
from ..util.colour_functions import complex_to_rgb
from ..util.constants import *


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

        if self.plot_type == 'volume':
            
            abs_max= np.amax(np.abs(eigenstates))
            psi = (psi)/(abs_max)

            L = self.eigenstates.extent/2/Å
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


        if self.plot_type == 'abs-volume':
            
            abs_max= np.amax(np.abs(eigenstates))
            psi = (psi)/(abs_max)

            L = self.eigenstates.extent/2/Å
            N = self.eigenstates.N

            vol = mlab.pipeline.volume(mlab.pipeline.scalar_field(np.abs(psi)), vmin= contrast_vals[0], vmax= contrast_vals[1])
            # Change the color transfer function

            mlab.outline()
            mlab.axes(xlabel='x [Å]', ylabel='y [Å]', zlabel='z [Å]',nb_labels=6 , ranges = (-L,L,-L,L,-L,L) )
            #azimuth angle
            φ = 30
            mlab.view(azimuth= φ,  distance=N*3.5)
            mlab.show()




        elif self.plot_type == 'contour':
            psi = eigenstates[k]
            L = self.eigenstates.extent/2/Å
            N = self.eigenstates.N
            isovalue = np.mean(contrast_vals)
            abs_max= np.amax(np.abs(eigenstates))
            psi = (psi)/(abs_max)

            field = mlab.pipeline.scalar_field(np.abs(psi))

            arr = mlab.screenshot(antialiased = False)

            mlab.outline()
            mlab.axes(xlabel='x [Å]', ylabel='y [Å]', zlabel='z [Å]',nb_labels=6 , ranges = (-L,L,-L,L,-L,L) )
            colour_data = np.angle(psi.T.ravel())%(2*np.pi)
            field.image_data.point_data.add_array(colour_data)
            field.image_data.point_data.get_array(1).name = 'phase'
            field.update()
            field2 = mlab.pipeline.set_active_attribute(field, 
                                                        point_scalars='scalar')
            contour = mlab.pipeline.contour(field2)
            contour.filter.contours= [isovalue,]
            contour2 = mlab.pipeline.set_active_attribute(contour, 
                                                        point_scalars='phase')
            s = mlab.pipeline.surface(contour, colormap='hsv', vmin= 0.0 ,vmax= 2.*np.pi)

            s.scene.light_manager.light_mode = 'vtk'
            s.actor.property.interpolation = 'phong'


            #azimuth angle
            φ = 30
            mlab.view(azimuth= φ,  distance=N*3.5)

            mlab.show()

    def animate(self,  contrast_vals= [0.1, 0.25]):
        eigenstates = self.eigenstates.array
        energies = self.eigenstates.energies
        mlab.figure(1, bgcolor=(0, 0, 0), size=(700, 700))

        
        if self.plot_type == 'volume':
            psi = eigenstates[0]
            
            abs_max= np.amax(np.abs(eigenstates))
            psi = (psi)/(abs_max)


            L = self.eigenstates.extent/2/Å
            N = self.eigenstates.N
            field = mlab.pipeline.scalar_field(psi)
            vol = mlab.pipeline.volume(field)

            color1 = complex_to_rgb(np.exp( 1j*2*np.pi/10*0)) 
            color2 = complex_to_rgb(-np.exp( 1j*2*np.pi/10*0)) 

            # Change the color transfer function
            from tvtk.util import ctf
            c = ctf.save_ctfs(vol._volume_property)
            c['rgb'] = [[-0.45, *color1],
                        [-0.4, *color1],
                        [-0.3, *color1],
                        [-0.2, *color1],
                        [-0.001, *color1],
                        [0.0, 0.0, 0.0, 0.0],
                        [0.001, *color2],
                        [0.2, *color2],
                        [0.3, *color2],
                        [0.4, *color2],
                        [0.45, *color2]]

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

                        color1 = complex_to_rgb(np.exp( 1j*2*np.pi/10*k1)*np.cos(np.pi*t) + np.exp( 1j*2*np.pi/10*k2)*np.sin(np.pi*t)) 
                        color2 = complex_to_rgb(-np.exp( 1j*2*np.pi/10*k1)*np.cos(np.pi*t) - np.exp( 1j*2*np.pi/10*k2)*np.sin(np.pi*t)) 
                    else:
                        psi = eigenstates[k1]
                        color1 = complex_to_rgb(np.exp( 1j*2*np.pi/10*k1)) 
                        color2 = complex_to_rgb(-np.exp( 1j*2*np.pi/10*k1)) 

                    psi = (psi)/(abs_max)
                    field.mlab_source.scalars = psi
                    # Change the color transfer function
                    from tvtk.util import ctf
                    c = ctf.save_ctfs(vol._volume_property)
                    c['rgb'] = [[-0.45, *color1],
                                [-0.4, *color1],
                                [-0.3, *color1],
                                [-0.2, *color1],
                                [-0.001, *color1],
                                [0.0, 0.0, 0.0, 0.0],
                                [0.001, *color2],
                                [0.2, *color2],
                                [0.3, *color2],
                                [0.4, *color2],
                                [0.45, *color2]]

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


        if self.plot_type == 'abs-volume':
            psi = eigenstates[0]
            
            abs_max= np.amax(np.abs(eigenstates))
            psi = np.abs((psi)/(abs_max))


            L = self.eigenstates.extent/2/Å
            N = self.eigenstates.N
            psi = np.where(psi > contrast_vals[1], contrast_vals[1],psi)
            psi = np.where(psi < contrast_vals[0], contrast_vals[0],psi)
            field = mlab.pipeline.scalar_field(psi)
            vol = mlab.pipeline.volume(field)


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

                    psi = np.abs((psi)/(abs_max))
                    psi = np.where(psi > contrast_vals[1], contrast_vals[1],psi)
                    psi = np.where(psi < contrast_vals[0], contrast_vals[0],psi)

                    field.mlab_source.scalars = psi
                    # Change the color transfer function

                    φ = 30 + data['t'] * 360 / 10 
                    mlab.view(azimuth= φ, distance=N*3.5)

                    yield

            animation()
            mlab.show()



        elif self.plot_type == 'contour':
            psi = eigenstates[0]
            L = self.eigenstates.extent/2/Å
            N = self.eigenstates.N
            isovalue = np.mean(contrast_vals)


            abs_max= np.amax(np.abs(eigenstates))
            psi = (psi)/(abs_max)

            field = mlab.pipeline.scalar_field(np.abs(psi))

            arr = mlab.screenshot(antialiased = False)

            mlab.outline()
            mlab.axes(xlabel='x [Å]', ylabel='y [Å]', zlabel='z [Å]',nb_labels=6 , ranges = (-L,L,-L,L,-L,L) )
            colour_data = np.angle(psi.T.ravel())%(2*np.pi)
            field.image_data.point_data.add_array(colour_data)
            field.image_data.point_data.get_array(1).name = 'phase'
            field.update()
            field2 = mlab.pipeline.set_active_attribute(field, 
                                                        point_scalars='scalar')
            contour = mlab.pipeline.contour(field2)
            contour.filter.contours= [isovalue,]
            contour2 = mlab.pipeline.set_active_attribute(contour, 
                                                        point_scalars='phase')
            s = mlab.pipeline.surface(contour2, colormap='hsv', vmin= 0.0 ,vmax= 2.*np.pi)

            s.scene.light_manager.light_mode = 'vtk'
            s.actor.property.interpolation = 'phong'


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
                        psi = (np.cos(np.pi*t)*eigenstates[k1]*np.exp( 1j*2*np.pi/10*k1) 
                             + np.sin(np.pi*t)*eigenstates[k2]*np.exp( 1j*2*np.pi/10*k2))


                    else:
                        psi = eigenstates[k1]*np.exp( 1j*2*np.pi/10*k1)
                    psi = (psi)/(abs_max)
                    np.copyto(colour_data, np.angle(psi.T.ravel())%(2*np.pi))
                    field.mlab_source.scalars = np.abs(psi)

                    φ = 30 + data['t'] * 360 / 10 
                    mlab.view(azimuth= φ, distance=N*3.5)


                    yield
            animation()
            mlab.show()







    def superpositions(self, states, contrast_vals= [0.1, 0.25], **kw):

        params = {'dt': 0.1}
        for k in kw.keys():
            if k in params:
                params[k] = kw[k]
            else:
                raise KeyError

        
        coeffs = states
        eigenstates = self.eigenstates.array
        energies = self.eigenstates.energies
        mlab.figure(1, bgcolor=(0, 0, 0), size=(700, 700))
        psi = sum([eigenstates[i]*coeffs[i] for i in range(len(coeffs))])

        if self.plot_type == 'volume':
            raise NotImplementedError
        elif self.plot_type == 'abs-volume':
            abs_max= np.amax(np.abs(eigenstates))
            psi = np.abs((psi)/(abs_max))


            L = self.eigenstates.extent/2/Å
            N = self.eigenstates.N
            psi = np.where(psi > contrast_vals[1], contrast_vals[1],psi)
            psi = np.where(psi < contrast_vals[0], contrast_vals[0],psi)
            field = mlab.pipeline.scalar_field(psi)
            vol = mlab.pipeline.volume(field)


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
                    data['t'] += params['dt']
                    t = data['t']
                    psi = sum([eigenstates[i]*np.exp(-1.0j*energies[i]*t)*coeffs[i]
                            for i in range(len(coeffs))])
                    psi = np.abs((psi)/(abs_max))

                    psi = np.where(psi > contrast_vals[1], contrast_vals[1],psi)
                    psi = np.where(psi < contrast_vals[0], contrast_vals[0],psi)
                    field.mlab_source.scalars = psi

                    φ = 30 + data['t'] * 360 / 10 
                    mlab.view(azimuth= φ, distance=N*3.5)

                    yield

            animation()
            mlab.show()
        elif self.plot_type == 'contour':
            L = self.eigenstates.extent/2/Å
            N = self.eigenstates.N
            isovalue = np.mean(contrast_vals)


            abs_max= np.amax(np.abs(eigenstates))
            psi = (psi)/(abs_max)

            field = mlab.pipeline.scalar_field(np.abs(psi))

            arr = mlab.screenshot(antialiased = False)

            mlab.outline()
            mlab.axes(xlabel='x [Å]', ylabel='y [Å]', zlabel='z [Å]',nb_labels=6 , ranges = (-L,L,-L,L,-L,L) )
            colour_data = np.angle(psi.T.ravel())%(2*np.pi)
            field.image_data.point_data.add_array(colour_data)
            field.image_data.point_data.get_array(1).name = 'phase'
            field.update()
            field2 = mlab.pipeline.set_active_attribute(field, 
                                                        point_scalars='scalar')
            contour = mlab.pipeline.contour(field2)
            contour.filter.contours= [isovalue,]
            contour2 = mlab.pipeline.set_active_attribute(contour, 
                                                        point_scalars='phase')
            s = mlab.pipeline.surface(contour2, colormap='hsv', vmin= 0.0 ,vmax= 2.*np.pi)

            s.scene.light_manager.light_mode = 'vtk'
            s.actor.property.interpolation = 'phong'


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
                               for i in range(len(coeffs))])

                    psi = (psi)/(abs_max)
                    np.copyto(colour_data, np.angle(psi.T.ravel())%(2*np.pi))
                    field.mlab_source.scalars = np.abs(psi)

                    φ = 30 + data['t'] * 360 / 10 
                    mlab.view(azimuth= φ, distance=N*3.5)
                    yield
            animation()
            mlab.show()

