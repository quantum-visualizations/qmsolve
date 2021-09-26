import numpy as np
import matplotlib.pyplot as plt
from matplotlib import widgets
from matplotlib import animation
from .visualization import Visualization
from ..util.constants import *

class VisualizationSingleParticle1D(Visualization):
    def __init__(self,eigenstates):
        self.eigenstates = eigenstates



    def plot_eigenstate(self, k, xlim = None, show_imaginary_part = False):
        eigenstates_array = self.eigenstates.array
        energies = self.eigenstates.energies
        plt.style.use("dark_background")

        fig = plt.figure(figsize=(16/9 *5.804 * 0.9,5.804)) 

        grid = plt.GridSpec(2, 2, width_ratios=[4.5, 1], height_ratios=[1, 1] , hspace=0.1, wspace=0.2)
        ax1 = fig.add_subplot(grid[0:2, 0:1])
        ax2 = fig.add_subplot(grid[0:2, 1:2])

        ax1.set_xlabel("x [Å]")

        ax2.set_title('Energy Level')
        ax2.set_facecolor('black')


        ax2.set_ylabel('$E_N$ [eV]')
        ax2.set_xticks(ticks=[])
        if xlim != None:
            ax1.set_xlim(np.array(xlim)/Å)

        ymax = np.amax(eigenstates_array)
        ymin = np.amin(eigenstates_array)
        ax1.set_ylim([ymin*1.2, ymax*1.2 ])



        E0 = energies[0]

        x = np.linspace(-self.eigenstates.extent/2, self.eigenstates.extent/2, self.eigenstates.N)


        if show_imaginary_part == True:
            eigenstate_plot1 = ax1.plot(x/Å, np.real(eigenstates_array[k]), label='$Re|\psi(x)|$')
            eigenstate_plot2 = ax1.plot(x/Å, np.imag(eigenstates_array[k]), label='$Im|\psi(x)|$')
            eigenstate_plot3 = ax1.plot(x/Å, np.abs(eigenstates_array[k]), label='$|\psi(x)|$', color='white')
            ax1.legend()
        else:
            eigenstate_plot1 = ax1.plot(x/Å, np.real(eigenstates_array[k]))


        for E in energies:
            ax2.plot([0,1], [E, E], color='gray', alpha=0.5)

        ax2.plot([0,1], [energies[k], energies[k]], color='yellow', lw = 3)
        plt.show()



    def slider_plot(self, xlim = None, show_imaginary_part = False):
        plt.style.use("dark_background")

        eigenstates_array = self.eigenstates.array
        energies = self.eigenstates.energies

        fig = plt.figure(figsize=(16/9 *5.804 * 0.9,5.804)) 

        grid = plt.GridSpec(2, 2, width_ratios=[5, 1], height_ratios=[1, 1] , hspace=0.1, wspace=0.2)
        ax1 = fig.add_subplot(grid[0:2, 0:1])
        ax2 = fig.add_subplot(grid[0:2, 1:2])

        ax1.set_xlabel("x [Å]")

        ax2.set_title('Energy Level')
        ax2.set_facecolor('black')

        ax2.set_ylabel('$E_N$ [eV]')
        ax2.set_xticks(ticks=[])
        if xlim != None:
            ax1.set_xlim(np.array(xlim)/Å)

        ymax = np.amax(eigenstates_array)
        ymin = np.amin(eigenstates_array)
        ax1.set_ylim([ymin*1.2, ymax*1.2 ])


        E0 = energies[0]
        for E in energies:
            ax2.plot([0,1], [E, E], color='gray', alpha=0.5)

        x = np.linspace(-self.eigenstates.extent/2, self.eigenstates.extent/2, self.eigenstates.N)

        if show_imaginary_part == True:
            eigenstate_plot1, = ax1.plot(x/Å, np.real(eigenstates_array[0]), label='$Re|\psi(x)|$')
            eigenstate_plot2, = ax1.plot(x/Å, np.imag(eigenstates_array[0]), label='$Im|\psi(x)|$')
            eigenstate_plot3, = ax1.plot(x/Å, np.abs(eigenstates_array[0]), label='$|\psi(x)|$', color='white')
            ax1.legend()
        else:
            eigenstate_plot, = ax1.plot(x/Å, np.real(eigenstates_array[0]))


        line = ax2.plot([0,1], [energies[1], energies[1]], color='yellow', lw = 3)

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
            if show_imaginary_part == True:

                eigenstate_plot1.set_ydata(np.real(eigenstates_array[state]))
                eigenstate_plot2.set_ydata(np.imag(eigenstates_array[state]))
                eigenstate_plot3.set_ydata( np.abs(eigenstates_array[state]))
            else:
                eigenstate_plot.set_ydata(np.real(eigenstates_array[state]))

            line[0].set_ydata([energies[state], energies[state]])

        slider.on_changed(update)
        plt.show()






    def animate(self,  seconds_per_eigenstate = 0.5, fps = 20, max_states = None, xlim = None, save_animation = False, show_imaginary_part = False):

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

        ax1.set_xlabel("[Å]")

        ax2.set_title('Energy Level')
        ax2.set_facecolor('black')

        ax2.set_ylabel('$E_N$ [eV]')
        ax2.set_xticks(ticks=[])
        if xlim != None:
            ax1.set_xlim(np.array(xlim)/Å)
        ymax = np.amax(eigenstates_array)
        ymin = np.amin(eigenstates_array)
        ax1.set_ylim([ymin*1.2, ymax*1.2 ])


        E0 = energies[0]
        for E in energies:
            ax2.plot([0,1], [E, E], color='gray', alpha=0.5)
        
        x = np.linspace(-self.eigenstates.extent/2, self.eigenstates.extent/2, self.eigenstates.N)
        if show_imaginary_part == True:
            eigenstate_plot1, = ax1.plot(x/Å, np.real(eigenstates_array[0]), label='$Re|\psi(x)|$')
            eigenstate_plot2, = ax1.plot(x/Å, np.imag(eigenstates_array[0]), label='$Im|\psi(x)|$')
            eigenstate_plot3, = ax1.plot(x/Å, np.abs(eigenstates_array[0]), label='$|\psi(x)|$', color='white')
            ax1.legend()
        else:
            eigenstate_plot, = ax1.plot(x/Å, np.real(eigenstates_array[0]))



        line, = ax2.plot([0,1], [energies[1], energies[1]], color='yellow', lw = 3)

        plt.subplots_adjust(bottom=0.2)

        import matplotlib.animation as animation

        animation_data = {'n': 0.0}
        def func_animation(*arg):
            animation_data['n'] = (animation_data['n'] + 0.1) % len(energies)
            state = int(animation_data['n'])
            if (animation_data['n'] % 1.0) > 0.5:
                transition_time = (animation_data['n'] - int(animation_data['n']) - 0.5)
                wavefunction = (np.cos(np.pi*transition_time)*eigenstates_array[state] + 
                                np.sin(np.pi*transition_time)*eigenstates_array[(state + 1) % len(energies)])

                if show_imaginary_part == True:

                    eigenstate_plot1.set_ydata(np.real(wavefunction))
                    eigenstate_plot2.set_ydata(np.imag(wavefunction))
                    eigenstate_plot3.set_ydata( np.abs(wavefunction))
                    return eigenstate_plot1, eigenstate_plot2, eigenstate_plot3, line
                else:
                    eigenstate_plot.set_ydata(np.real(wavefunction))
                    return eigenstate_plot, line

                E_N = energies[state] 
                E_M = energies[(state + 1) % len(energies)]
                E =  E_N*np.cos(np.pi*transition_time)**2 + E_M*np.sin(np.pi*transition_time)**2
                line.set_ydata([E, E])
            else:
                line.set_ydata([energies[state], energies[state]])
                wavefunction = eigenstates_array[int(state)]

                if show_imaginary_part == True:

                    eigenstate_plot1.set_ydata(np.real(wavefunction))
                    eigenstate_plot2.set_ydata(np.imag(wavefunction))
                    eigenstate_plot3.set_ydata( np.abs(wavefunction))
                    return eigenstate_plot1, eigenstate_plot2, eigenstate_plot3, line
                else:
                    eigenstate_plot.set_ydata(np.real(wavefunction))
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
        """
        Visualize the time evolution of a superposition of energy eigenstates.
        The circle widgets control the relative phase of each of the eigenstates.
        These widgets are inspired by the circular phasors from the
        quantum mechanics applets by Paul Falstad:
        https://www.falstad.com/qm1d/
        """

        total_frames = fps * total_time
        from .complex_slider_widget import ComplexSliderWidget
        eigenstates = self.eigenstates.array
        coeffs = None
        get_norm_factor = lambda psi: 1.0/np.sqrt(np.sum(psi*np.conj(psi)))
        animation_data = {'ticks': 0, 'norm': get_norm_factor(eigenstates[0]),
                          'is_paused': False}
        psi0 = eigenstates[0]*get_norm_factor(eigenstates[0])
        if isinstance(states, int) or isinstance(states, float):
            coeffs = np.array([1.0 if i == 0 else 0.0 for i in range(states)],
                           dtype=np.complex128)
            eigenstates = eigenstates[0: states]
        else:
            coeffs = states
            eigenstates = eigenstates[0: len(states)]
            states = len(states)
            psi0 = np.tensordot(coeffs, eigenstates, 1)
            animation_data['norm'] = get_norm_factor(psi0)
            psi0 *= animation_data['norm']
        energies = self.eigenstates.energies
        params = {'dt': 0.001, 
                  'xlim': [-self.eigenstates.extent/2.0, 
                         self.eigenstates.extent/2.0],
                  'save_animation': False,
                  'frames': 120
                 }
        for k in kw.keys():
            params[k] = kw[k]

        plt.style.use("dark_background")
        fig = plt.figure(figsize=(16/9 *5.804 * 0.9,5.804)) 
        grid = plt.GridSpec(5, states)
        ax = fig.add_subplot(grid[0:3, 0:states])
        
        ax.set_xlabel("[Å]")
        x = np.linspace(-self.eigenstates.extent/2.0,
                        self.eigenstates.extent/2.0,
                        len(eigenstates[0]))
        ax.set_yticks([])
        ax.set_xlim(np.array(params['xlim'])/Å)

        line1, = ax.plot(x/Å, np.real(eigenstates[0]), label='$Re|\psi(x)|$')
        line2, = ax.plot(x/Å, np.imag(eigenstates[0]), label='$Im|\psi(x)|$')
        line3, = ax.plot(x/Å, np.abs(eigenstates[0]), label='$|\psi(x)|$', color='white')
        ax.set_ylim(-1.7*np.amax(np.abs(psi0)), 1.7*np.amax(np.abs(psi0)))
        ax.legend()

        def make_update(n):
            def update(phi, r):
                animation_data['is_paused'] = True
                coeffs[n] = r*np.exp(1.0j*phi)
                psi = np.tensordot(coeffs, eigenstates, 1)
                animation_data['norm'] = get_norm_factor(psi)
                line1.set_ydata(np.real(psi))
                line2.set_ydata(np.imag(psi))
                line3.set_ydata(np.abs(psi))
            return update

        widgets = []
        circle_artists = []
        for i in range(states):
            circle_ax = fig.add_subplot(grid[4, i], projection='polar')
            circle_ax.set_title('n=' + str(i) # + '\nE=' + str() + '$E_0$'
                                )
            circle_ax.set_xticks([])
            circle_ax.set_yticks([])
            widgets.append(ComplexSliderWidget(circle_ax, 0.0, 1.0, animated=True))
            widgets[i].on_changed(make_update(i))
            circle_artists.append(widgets[i].get_artist())
        artists = circle_artists + [line1, line2, line3]

        def func(*args):
            animation_data['ticks'] += 1
            e = 1.0
            if animation_data['is_paused']:
                animation_data['is_paused'] = False
            else:
                e *= np.exp(-1.0j*energies[0:states]*params['dt'])
            np.copyto(coeffs, coeffs*e)
            norm_factor = animation_data['norm']
            psi = np.tensordot(coeffs*norm_factor, eigenstates, 1)
            line1.set_ydata(np.real(psi))
            line2.set_ydata(np.imag(psi))
            line3.set_ydata(np.abs(psi))
            if animation_data['ticks'] % 2:
                return [line1, line2, line3]
            else:
                for i, c in enumerate(coeffs):
                    phi, r = np.angle(c), np.abs(c)
                    artists[i].set_xdata([phi, phi])
                    artists[i].set_ydata([0.0, r])
                return artists
        a = animation.FuncAnimation(fig, func, blit=True, interval=1000.0/60.0,
                                    frames=None if (not params['save_animation']) else
                                    total_frames)
        if params['save_animation'] == True:
            Writer = animation.writers['ffmpeg']
            writer = Writer(fps=fps, metadata=dict(artist='Me'), 
                            bitrate=-1)
            a.save('animation.mp4', writer=writer)
            return
        plt.show()





from .visualization import TimeVisualization

class TimeVisualizationSingleParticle1D(TimeVisualization):
    def __init__(self,simulation):
        self.simulation = simulation
        self.H = simulation.H

    def plot(self, t, xlim=None, figsize=(16/9 *5.804 * 0.9,5.804)):


        self.simulation.Ψ_plot = self.simulation.Ψ/self.simulation.Ψmax
        plt.style.use("dark_background")

        fig = plt.figure(figsize=figsize)
        
        ax = fig.add_subplot(1, 1, 1)

        ax.set_xlabel("[Å]")
        ax.set_title("$\psi(x,t)$")

        time_ax = ax.text(0.97,0.97, "",  color = "white",
                        transform=ax.transAxes, ha="right", va="top")
        time_ax.set_text(u"t = {} femtoseconds".format("%.3f"  % (t/femtoseconds)))



        if xlim != None:
            ax.set_xlim(np.array(xlim)/Å)


        index = int((self.simulation.store_steps)/self.simulation.total_time*t)
        
        L = self.simulation.H.extent/Å

        x = np.linspace(-self.simulation.H.extent/2, self.simulation.H.extent/2, self.simulation.H.N)


        potential_plot = ax.plot(x/Å, (self.simulation.H.Vgrid + self.simulation.Vmin)/(self.simulation.Vmax-self.simulation.Vmin), label='$V(x)$')  
        real_plot = ax.plot(x/Å, np.real(self.simulation.Ψ_plot[index]), label='$Re|\psi(x)|$')
        imag_plot = ax.plot(x/Å, np.imag(self.simulation.Ψ_plot[index]), label='$Im|\psi(x)|$')
        abs_plot = ax.plot(x/Å, np.abs(self.simulation.Ψ_plot[index]), label='$|\psi(x)|$', color='white')
        ax.legend('lower left')

        plt.show()


    def animate(self,xlim=None, figsize=(16/9 *5.804 * 0.9,5.804), animation_duration = 5, fps = 20, save_animation = False):
        total_frames = int(fps * animation_duration)
        dt = self.simulation.total_time/total_frames
        self.simulation.Ψ_plot = self.simulation.Ψ/self.simulation.Ψmax
        plt.style.use("dark_background")

        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(1, 1, 1)
        index = 0
        
        ax.set_xlabel("[Å]")
        ax.set_title("$\psi(x,t)$")

        time_ax = ax.text(0.97,0.97, "",  color = "white",
                        transform=ax.transAxes, ha="right", va="top")
        time_ax.set_text(u"t = {} femtoseconds".format("%.3f"  % (0./femtoseconds)))



        if xlim != None:
            ax.set_xlim(np.array(xlim)/Å)


        index = 0
        
        L = self.simulation.H.extent/Å

        x = np.linspace(-self.simulation.H.extent/2, self.simulation.H.extent/2, self.simulation.H.N)


        potential_plot = ax.plot(x/Å, (self.simulation.H.Vgrid + self.simulation.Vmin)/(self.simulation.Vmax-self.simulation.Vmin), label='$V(x)$')  
        real_plot, = ax.plot(x/Å, np.real(self.simulation.Ψ_plot[index]), label='$Re|\psi(x)|$')
        imag_plot, = ax.plot(x/Å, np.imag(self.simulation.Ψ_plot[index]), label='$Im|\psi(x)|$')
        abs_plot, = ax.plot(x/Å, np.abs(self.simulation.Ψ_plot[index]), label='$|\psi(x)|$', color='white')
        ax.legend(loc = 'lower left')

        import matplotlib.animation as animation


        #print(total_frames)
        animation_data = {'t': 0.0, 'x' :x, 'ax':ax ,'frame' : 0}
        def func_animation(*arg):
            
            time_ax.set_text(u"t = {} femtoseconds".format("%.3f"  % (animation_data['t']/femtoseconds)))

            animation_data['t'] = animation_data['t'] + dt
            if animation_data['t'] > self.simulation.total_time:
                animation_data['t'] = 0.0

            #print(animation_data['frame'])
            animation_data['frame'] +=1
            index = int((self.simulation.store_steps)/self.simulation.total_time * animation_data['t'])

            real_plot.set_ydata(np.real(self.simulation.Ψ_plot[index]))
            imag_plot.set_ydata(np.imag(self.simulation.Ψ_plot[index]))
            abs_plot.set_ydata(np.abs(self.simulation.Ψ_plot[index]))


            return potential_plot,real_plot,imag_plot,abs_plot, time_ax

        frame = 0
        a = animation.FuncAnimation(fig, func_animation,
                                    blit=False, frames=total_frames, interval= 1/fps * 1000)
        if save_animation == True:
            Writer = animation.writers['ffmpeg']
            writer = Writer(fps=fps, metadata=dict(artist='Me'), bitrate=1800)
            a.save('animation.mp4', writer=writer)
        else:
            plt.show()
