import numpy as np
import matplotlib.pyplot as plt
from matplotlib import widgets
from matplotlib import animation
from .visualization import Visualization
from ..util.colour_functions import complex_to_rgb, complex_to_rgba
from ..util.constants import *


class VisualizationSingleParticle2D(Visualization):
    def __init__(self,eigenstates):
        self.eigenstates = eigenstates



    def plot_eigenstate(self, k, xlim=None, ylim=None):
        eigenstates_array = self.eigenstates.array
        energies = self.eigenstates.energies

        plt.style.use("dark_background")
        fig = plt.figure(figsize=(16/9 *5.804 * 0.9,5.804)) 

        grid = plt.GridSpec(2, 2, width_ratios=[4.5, 1], height_ratios=[1, 1] , hspace=0.1, wspace=0.2)
        ax1 = fig.add_subplot(grid[0:2, 0:1])
        ax2 = fig.add_subplot(grid[0:2, 1:2])

        ax1.set_xlabel("$x$ [Å]")
        ax1.set_ylabel("$y$ [Å]")
        ax1.set_title("$\Psi(x,y)$")

        ax2.set_title('Energy Level')
        ax2.set_facecolor('black')
        ax2.set_ylabel('$E_N$ [eV]')
        ax2.set_xticks(ticks=[])

        if xlim != None:
            ax1.set_xlim(np.array(xlim)/Å)
        if ylim != None:
            ax1.set_ylim(np.array(ylim)/Å)

        E0 = energies[0]

        for E in energies:
            ax2.plot([0,1], [E, E], color='gray', alpha=0.5)

        ax2.plot([0,1], [energies[k], energies[k]], color='yellow', lw = 3)

        ax1.set_aspect('equal')
        L =  self.eigenstates.extent/2/Å
        im = ax1.imshow(complex_to_rgb(eigenstates_array[k]*np.exp( 1j*2*np.pi/10*k)), origin='lower',extent = [-L, L, -L, L],  interpolation = 'bilinear')
        plt.show()


    def slider_plot(self, xlim=None, ylim=None):

        eigenstates_array = self.eigenstates.array
        energies = self.eigenstates.energies

        plt.style.use("dark_background")
        fig = plt.figure(figsize=(16/9 *5.804 * 0.9,5.804)) 

        grid = plt.GridSpec(2, 2, width_ratios=[5, 1], height_ratios=[1, 1] , hspace=0.1, wspace=0.2)
        ax1 = fig.add_subplot(grid[0:2, 0:1])
        ax2 = fig.add_subplot(grid[0:2, 1:2])

        ax1.set_xlabel("$x$ [Å]")
        ax1.set_ylabel("$y$ [Å]")
        ax1.set_title("$\Psi(x,y)$")

        ax2.set_title('Energy Level')
        ax2.set_facecolor('black')
        ax2.set_ylabel('$E_N$ [eV]')
        ax2.set_xticks(ticks=[])

        if xlim != None:
            ax1.set_xlim(np.array(xlim)/Å)
        if ylim != None:
            ax1.set_ylim(np.array(ylim)/Å)


        E0 = energies[0]
        for E in energies:
            ax2.plot([0,1], [E, E], color='gray', alpha=0.5)


        ax1.set_aspect('equal')
        L = self.eigenstates.extent/2/Å
        eigenstate_plot = ax1.imshow(complex_to_rgb(eigenstates_array[0]*np.exp( 1j*2*np.pi/10*0)), origin='lower',extent = [-L, L, -L, L],  interpolation = 'bilinear')
        
        line = ax2.plot([0,1], [energies[0], energies[0]], color='yellow', lw = 3)

        plt.subplots_adjust(bottom=0.2)
        from matplotlib.widgets import Slider
        slider_ax = plt.axes([0.2, 0.05, 0.7, 0.05])
        slider = Slider(slider_ax,      # the axes object containing the slider
                          'state',            # the name of the slider parameter
                          0,          # minimal value of the parameter
                          len(eigenstates_array)-1,          # maximal value of the parameter
                          valinit = 0,  # initial value of the parameter 
                          valstep = 1,
                          color = '#5c05ff' 
                         )

        def update(state):
            state = int(state)
            eigenstate_plot.set_data(complex_to_rgb(eigenstates_array[state]*np.exp( 1j*2*np.pi/10*state)))
            line[0].set_ydata([energies[state], energies[state]])

        slider.on_changed(update)
        plt.show()






    def animate(self,seconds_per_eigenstate = 0.5, fps = 20, max_states = None, xlim=None, ylim=None, save_animation = False):

        if max_states == None:
            max_states = len(self.eigenstates.energies)

        frames_per_eigenstate = fps * seconds_per_eigenstate
        total_time = max_states * seconds_per_eigenstate
        total_frames = int(fps * total_time)


        eigenstates_array = self.eigenstates.array
        energies = self.eigenstates.energies

        plt.style.use("dark_background")
        fig = plt.figure(figsize=(16/9 *5.804 * 0.9,5.804)) 

        grid = plt.GridSpec(2, 2, width_ratios=[5, 1], height_ratios=[1, 1] , hspace=0.1, wspace=0.2)
        ax1 = fig.add_subplot(grid[0:2, 0:1])
        ax2 = fig.add_subplot(grid[0:2, 1:2])

        ax1.set_xlabel("$x$ [Å]")
        ax1.set_ylabel("$y$ [Å]")
        ax1.set_title("$\Psi(x,y)$")

        ax2.set_title('Energy Level')
        ax2.set_facecolor('black')
        ax2.set_ylabel('$E_N$ [eV]')
        ax2.set_xticks(ticks=[])

        if xlim != None:
            ax1.set_xlim(np.array(xlim)/Å)
        if ylim != None:
            ax1.set_ylim(np.array(ylim)/Å)

        E0 = energies[0]
        for E in energies:
            ax2.plot([0,1], [E, E], color='gray', alpha=0.5)
        
        ax1.set_aspect('equal')
        L = self.eigenstates.extent/2/Å
        eigenstate_plot = ax1.imshow(complex_to_rgb(eigenstates_array[0]*np.exp( 1j*2*np.pi/10*0)),  origin='lower',extent = [-L, L, -L, L],   interpolation = 'bilinear')

        line, = ax2.plot([0,1], [energies[0], energies[0]], color='yellow', lw = 3)

        plt.subplots_adjust(bottom=0.2)

        import matplotlib.animation as animation

        animation_data = {'n': 0.0}
        Δn = 1/frames_per_eigenstate

        def func_animation(*arg):
            animation_data['n'] = (animation_data['n'] + Δn) % len(energies)
            state = int(animation_data['n'])
            if (animation_data['n'] % 1.0) > 0.5:
                transition_time = (animation_data['n'] - int(animation_data['n']) - 0.5)
                eigenstate_combination = (np.cos(np.pi*transition_time)*eigenstates_array[state]*np.exp( 1j*2*np.pi/10*state) + 
                                         np.sin(np.pi*transition_time)*
                                         eigenstates_array[(state + 1) % len(energies)]*np.exp( 1j*2*np.pi/10*(state + 1)) )
                
                eigenstate_plot.set_data(complex_to_rgb(eigenstate_combination))


                E_N = energies[state] 
                E_M = energies[(state + 1) % len(energies)]
                E =  E_N*np.cos(np.pi*transition_time)**2 + E_M*np.sin(np.pi*transition_time)**2
                line.set_ydata([E, E])
            else:
                line.set_ydata([energies[state], energies[state]])
                eigenstate_combination = eigenstates_array[int(state)]*np.exp( 1j*2*np.pi/10*state)
                eigenstate_plot.set_data(complex_to_rgb(eigenstate_combination))
            return eigenstate_plot, line

        a = animation.FuncAnimation(fig, func_animation,
                                    blit=True, frames=total_frames, interval= 1/fps * 1000)
        if save_animation == True:
            Writer = animation.writers['ffmpeg']
            writer = Writer(fps=fps, metadata=dict(artist='Me'), bitrate=1800)
            a.save('animation.mp4', writer=writer)
        else:
            plt.show()


    def superpositions(self, states, fps = 30, total_time = 20, **kw):
        params = {'dt': 0.001, 'xlim': [-self.eigenstates.extent/2, 
                                        self.eigenstates.extent/2],
                  'ylim': [-self.eigenstates.extent/2, 
                           self.eigenstates.extent/2],
                  'save_animation': False,
                  'hide_controls': False,
                  # 'plot_style': 'dark_background'
                 }
        for k in kw.keys():
            params[k] = kw[k]
        total_frames = fps * total_time
        from .complex_slider_widget import ComplexSliderWidget
        eigenstates = self.eigenstates.array
        energies = self.eigenstates.energies
        eigenstates = np.array(eigenstates)
        energies = np.array(energies)
        coeffs = None
        if isinstance(states, int) or isinstance(states, float):
            coeffs = np.array([1.0 if i == 0 else 0.0 for i in range(states)],
                           dtype=np.complex128)
            eigenstates = eigenstates[0: states]
        else:
            coeffs = states
            eigenstates = eigenstates[0: len(states)]
            states = len(states)
            params[k] = kw[k]
        N = eigenstates.shape[1]
        plt.style.use("dark_background")
        fig = plt.figure(figsize=(16/9 *5.804 * 0.9,5.804))
        grid_width = 10
        grid_length = states if states < 30 else 30
        grid = plt.GridSpec(grid_width, grid_length)
        grid_slice = grid[0:int(0.7*grid_width), 0:grid_length]
        if params['hide_controls']:
            grid_slice = grid[0:grid_width, 0:grid_length]
        ax = fig.add_subplot(grid_slice)
        ax.set_title("$\psi(x, y)$")
        ax.set_xlabel("$x$ [Å]")
        ax.set_ylabel("$y$ [Å]")
        # ax.set_xticks([])
        # ax.set_yticks([])
        get_norm_factor = lambda psi: 1.0/np.sqrt(np.sum(psi*np.conj(psi)))
        coeffs = np.array(coeffs, dtype=np.complex128)
        X, Y = np.meshgrid(np.linspace(-1.0, 1.0, eigenstates[0].shape[0]),
                        np.linspace(-1.0, 1.0, eigenstates[0].shape[1]))
        maxval = np.amax(np.abs(eigenstates[0]))


        ax.set_xlim(np.array(params['xlim'])/Å)
        ax.set_ylim(np.array(params['ylim'])/Å)


        im = plt.imshow(complex_to_rgb(eigenstates[0]), interpolation='bilinear',
                        origin='lower', extent=[-self.eigenstates.extent/2/Å, 
                                                self.eigenstates.extent/2/Å,
                                                -self.eigenstates.extent/2/Å, 
                                                self.eigenstates.extent/2/Å]
                        )
        # im2 = plt.imshow(0.0*eigenstates[0], cmap='gray')
        animation_data = {'ticks': 0, 'norm': 1.0}

        def make_update(n):
            def update(phi, r):
                coeffs[n] = r*np.exp(1.0j*phi)
                psi = np.dot(coeffs, 
                             eigenstates.reshape([states, N*N]))
                psi = psi.reshape([N, N])
                animation_data['norm'] = get_norm_factor(psi)
                psi *= animation_data['norm']
                # apsi = np.abs(psi)
                # im.set_alpha(apsi/np.amax(apsi))
            return update

        widgets = []
        circle_artists = []
        if not params['hide_controls']:
            for i in range(states):
                if states <= 30:
                    circle_ax = fig.add_subplot(grid[8:10, i], projection='polar')
                    circle_ax.set_title('n=' + str(i) # + '\nE=' + str() + '$E_0$'
                                        , size=8.0 if states < 15 else 6.0 
                                        )
                else:
                    circle_ax = fig.add_subplot(grid[8 if i < 30 else 9,
                                                     i if i < 30 else i-30], 
                                                projection='polar')
                    circle_ax.set_title('n=' + str(i) # + '\nE=' + str() + '$E_0$'
                                        , size=8.0 if states < 15 else 6.0 
                                        )
                circle_ax.set_xticks([])
                circle_ax.set_yticks([])
                widgets.append(ComplexSliderWidget(circle_ax, 0.0, 1.0, animated=True))
                widgets[i].on_changed(make_update(i))
                circle_artists.append(widgets[i].get_artist())
        artists = circle_artists + [im]

        def func(*args):
            animation_data['ticks'] += 1
            e = np.exp(-1.0j*energies[0:states]*params['dt'])
            np.copyto(coeffs, coeffs*e)
            norm_factor = animation_data['norm']
            psi = np.dot(coeffs*norm_factor, 
                         eigenstates.reshape([
                            states, N*N]))
            psi = psi.reshape([N, N])
            im.set_data(complex_to_rgb(psi))
            # apsi = np.abs(psi)
            # im.set_alpha(apsi/np.amax(apsi))
            # if animation_data['ticks'] % 2:
            #     return (im, )
            # else:
            if not params['hide_controls']:
                for i, c in enumerate(coeffs):
                    phi, r = np.angle(c), np.abs(c)
                    artists[i].set_xdata([phi, phi])
                    artists[i].set_ydata([0.0, r])
            return artists

        a = animation.FuncAnimation(fig, func, blit=True, interval= 1/fps * 1000,
                                    frames=None if (not params['save_animation']) else
                                    total_frames)
        if params['save_animation'] == True:
            Writer = animation.writers['ffmpeg']
            writer = Writer(fps=fps, metadata=dict(artist='Me'), 
                            bitrate=-1)
            a.save('animation.mp4', writer=writer)
            return
        plt.show()
        plt.show()



from .visualization import TimeVisualization

class TimeVisualizationSingleParticle2D(TimeVisualization):
    def __init__(self,simulation):
        self.simulation = simulation
        self.H = simulation.H

    def plot(self, t, xlim=None, ylim=None, figsize=(7, 7), potential_saturation=0.8, wavefunction_saturation=1.0):


        self.simulation.Ψ_plot = self.simulation.Ψ/self.simulation.Ψmax
        plt.style.use("dark_background")

        fig = plt.figure(figsize=figsize)
        
        ax = fig.add_subplot(1, 1, 1)

        ax.set_xlabel("[Å]")
        ax.set_ylabel("[Å]")
        ax.set_title("$\psi(x,y,t)$")

        time_ax = ax.text(0.97,0.97, "",  color = "white",
                        transform=ax.transAxes, ha="right", va="top")
        time_ax.set_text(u"t = {} femtoseconds".format("%.3f"  % (t/femtoseconds)))



        if xlim != None:
            ax.set_xlim(np.array(xlim)/Å)
        if ylim != None:
            ax.set_ylim(np.array(ylim)/Å)


        from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
        index = int((self.simulation.store_steps)/self.simulation.total_time*t)
        
        L = self.simulation.H.extent/Å
        ax.imshow((self.simulation.H.Vgrid + self.simulation.Vmin)/(self.simulation.Vmax-self.simulation.Vmin), vmax = 1.0/potential_saturation, vmin = 0, cmap = "gray", origin = "lower", interpolation = "bilinear", extent = [-L/2, L/2, -L/2, L/2])  

        ax.imshow(complex_to_rgba(self.simulation.Ψ_plot[index], max_val= wavefunction_saturation), origin = "lower", interpolation = "bilinear", extent = [-L/2, L/2, -L/2, L/2])  
        plt.show()


    def animate(self,xlim=None, ylim=None, figsize=(7, 7), animation_duration = 5, fps = 20, save_animation = False, potential_saturation=0.8, wavefunction_saturation=0.8):
        total_frames = int(fps * animation_duration)
        dt = self.simulation.total_time/total_frames
        self.simulation.Ψ_plot = self.simulation.Ψ/self.simulation.Ψmax
        plt.style.use("dark_background")

        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(1, 1, 1)
        index = 0
        
        L = self.simulation.H.extent/Å
        potential_plot = ax.imshow((self.simulation.H.Vgrid + self.simulation.Vmin)/(self.simulation.Vmax-self.simulation.Vmin), vmax = 1.0/potential_saturation, vmin = 0, cmap = "gray", origin = "lower", interpolation = "bilinear", extent = [-L/2, L/2, -L/2, L/2])  
        wavefunction_plot = ax.imshow(complex_to_rgba(self.simulation.Ψ_plot[0], max_val= wavefunction_saturation), origin = "lower", interpolation = "bilinear", extent=[-L / 2,L / 2,-L / 2,L / 2])


        if xlim != None:
            ax.set_xlim(np.array(xlim)/Å)
        if ylim != None:
            ax.set_ylim(np.array(ylim)/Å)

        import matplotlib.animation as animation

        ax.set_title("$\psi(x,y,t)$")
        ax.set_xlabel('[Å]')
        ax.set_ylabel('[Å]')

        time_ax = ax.text(0.97,0.97, "",  color = "white",
                        transform=ax.transAxes, ha="right", va="top")
        time_ax.set_text(u"t = {} femtoseconds".format("%.3f"  % 0.00))


        #print(total_frames)
        animation_data = {'t': 0.0, 'ax':ax ,'frame' : 0}
        def func_animation(*arg):
            
            time_ax.set_text(u"t = {} femtoseconds".format("%.3f"  % (animation_data['t']/femtoseconds)))

            animation_data['t'] = animation_data['t'] + dt
            if animation_data['t'] > self.simulation.total_time:
                animation_data['t'] = 0.0

            #print(animation_data['frame'])
            animation_data['frame'] +=1
            index = int((self.simulation.store_steps)/self.simulation.total_time * animation_data['t'])

            wavefunction_plot.set_data(complex_to_rgba(self.simulation.Ψ_plot[index], max_val= wavefunction_saturation))
            return potential_plot,wavefunction_plot, time_ax

        frame = 0
        a = animation.FuncAnimation(fig, func_animation,
                                    blit=True, frames=total_frames, interval= 1/fps * 1000)
        if save_animation == True:
            Writer = animation.writers['ffmpeg']
            writer = Writer(fps=fps, metadata=dict(artist='Me'), bitrate=1800)
            a.save('animation.mp4', writer=writer)
        else:
            plt.show()
