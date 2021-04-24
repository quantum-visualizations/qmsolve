import numpy as np
import matplotlib.pyplot as plt
from matplotlib import widgets
from matplotlib import animation
from .visualization import Visualization
from ..util.colour_functions import complex_to_rgb


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

        ax2.set_title('E Level')
        ax2.set_facecolor('black')
        ax2.set_ylabel('$E_N$ (Relative to $E_{1}$)')
        ax2.set_xticks(ticks=[])

        if xlim != None:
            ax1.set_xlim(xlim)
        if ylim != None:
            ax1.set_ylim(ylim)

        E0 = energies[0]

        for E in energies:
            ax2.plot([0,1], [E/E0, E/E0], color='gray', alpha=0.5)

        ax2.plot([0,1], [energies[k]/E0, energies[k]/E0], color='yellow', lw = 3)

        ax1.set_aspect('equal')
        L =  self.eigenstates.extent/2
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

        ax2.set_title('E Level')
        ax2.set_facecolor('black')
        ax2.set_ylabel('$E_N$ (Relative to $E_{1}$)')
        ax2.set_xticks(ticks=[])

        if xlim != None:
            ax1.set_xlim(xlim)
        if ylim != None:
            ax1.set_ylim(ylim)


        E0 = energies[0]
        for E in energies:
            ax2.plot([0,1], [E/E0, E/E0], color='gray', alpha=0.5)


        ax1.set_aspect('equal')
        L = self.eigenstates.extent/2
        eigenstate_plot = ax1.imshow(complex_to_rgb(eigenstates_array[0]*np.exp( 1j*2*np.pi/10*0)), origin='lower',extent = [-L, L, -L, L],  interpolation = 'bilinear')
        
        line = ax2.plot([0,1], [energies[0]/E0, energies[0]/E0], color='yellow', lw = 3)

        plt.subplots_adjust(bottom=0.2)
        from matplotlib.widgets import Slider
        slider_ax = plt.axes([0.2, 0.05, 0.7, 0.05])
        slider = Slider(slider_ax,      # the axes object containing the slider
                          'state',            # the name of the slider parameter
                          0,          # minimal value of the parameter
                          len(eigenstates_array)-1,          # maximal value of the parameter
                          valinit = 1,  # initial value of the parameter 
                          valstep = 1,
                          color = '#5c05ff' 
                         )

        def update(state):
            state = int(state)
            eigenstate_plot.set_data(complex_to_rgb(eigenstates_array[state]*np.exp( 1j*2*np.pi/10*state)))
            line[0].set_ydata([energies[state]/E0, energies[state]/E0])

        slider.on_changed(update)
        plt.show()






    def animate(self, xlim=None, ylim=None):

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

        ax2.set_title('E Level')
        ax2.set_facecolor('black')
        ax2.set_ylabel('$E_N$ (Relative to $E_{1}$)')
        ax2.set_xticks(ticks=[])

        if xlim != None:
            ax1.set_xlim(xlim)
        if ylim != None:
            ax1.set_ylim(ylim)

        E0 = energies[0]
        for E in energies:
            ax2.plot([0,1], [E/E0, E/E0], color='gray', alpha=0.5)
        
        # ax1.set_xlim( ??? )
        ax1.set_aspect('equal')
        L = self.eigenstates.extent/2
        eigenstate_plot = ax1.imshow(complex_to_rgb(eigenstates_array[0]*np.exp( 1j*2*np.pi/10*0)),  origin='lower',extent = [-L, L, -L, L],   interpolation = 'bilinear')

        line, = ax2.plot([0,1], [energies[0]/E0, energies[0]/E0], color='yellow', lw = 3)

        plt.subplots_adjust(bottom=0.2)

        import matplotlib.animation as animation

        animation_data = {'n': 0.0}
        def func_animation(*arg):
            animation_data['n'] = (animation_data['n'] + 0.1) % len(energies)
            state = int(animation_data['n'])
            if (animation_data['n'] % 1.0) > 0.5:
                transition_time = (animation_data['n'] - int(animation_data['n']) - 0.5)
                eigenstate_combination = (np.cos(np.pi*transition_time)*eigenstates_array[state]*np.exp( 1j*2*np.pi/10*state) + 
                                         np.sin(np.pi*transition_time)*
                                         eigenstates_array[(state + 1) % len(energies)]*np.exp( 1j*2*np.pi/10*(state + 1)) )
                
                eigenstate_plot.set_data(complex_to_rgb(eigenstate_combination))


                E_N = energies[state]/E0 
                E_M = energies[(state + 1) % len(energies)]/E0
                E =  E_N*np.cos(np.pi*transition_time)**2 + E_M*np.sin(np.pi*transition_time)**2
                line.set_ydata([E, E])
            else:
                line.set_ydata([energies[state]/E0, energies[state]/E0])
                eigenstate_combination = eigenstates_array[int(state)]*np.exp( 1j*2*np.pi/10*state)
                eigenstate_plot.set_data(complex_to_rgb(eigenstate_combination))
            return eigenstate_plot, line

        a = animation.FuncAnimation(fig, func_animation,
                                    blit=True, interval=1.0)
        plt.show()
        """
        # save animation
        Writer = animation.writers['ffmpeg']
        writer = Writer(fps=20, metadata=dict(artist='Me'), bitrate=1800)
        a.save('im.gif', writer=writer)
        """

    def superpositions(self, states, **kw):
        from .complex_slider_widget import ComplexSliderWidget
        eigenstates = self.eigenstates.array
        energies = self.eigenstates.energies
        eigenstates = np.array(eigenstates)
        energies = np.array(energies)
        params = {'dt': 0.001}
        for k in kw.keys():
            params[k] = kw[k]
        N = eigenstates.shape[1]

        plt.style.use("dark_background")
        fig = plt.figure(figsize=(16/9 *5.804 * 0.9,5.804)) 
        grid = plt.GridSpec(4, 10)
        ax = fig.add_subplot(grid[0:3, 0:10])
        ax.set_xticks([])
        ax.set_yticks([])
        get_norm_factor = lambda psi: 1.0/np.sqrt(np.sum(psi*np.conj(psi)))
        coeffs = np.array([1.0 if i == 0 else 0.0 for i in range(states)],
                        dtype=np.complex128)
        X, Y = np.meshgrid(np.linspace(-1.0, 1.0, eigenstates[0].shape[0]),
                        np.linspace(-1.0, 1.0, eigenstates[0].shape[1]))
        maxval = np.amax(np.abs(eigenstates[0]))
        im = plt.imshow(complex_to_rgb(eigenstates[0]), interpolation='bilinear')
        im2 = plt.imshow(0.0*eigenstates[0], cmap='gray')
        animation_data = {'ticks': 0, 'norm': 1.0}

        def make_update(n):
            def update(phi, r):
                coeffs[n] = r*np.exp(1.0j*phi)
                psi = np.dot(coeffs, 
                            eigenstates[0:states].reshape([states, N*N]))
                psi = psi.reshape([N, N])
                animation_data['norm'] = get_norm_factor(psi)
                psi *= animation_data['norm']
                # apsi = np.abs(psi)
                # im.set_alpha(apsi/np.amax(apsi))
            return update

        widgets = []
        circle_artists = []
        for i in range(states):
            circle_ax = fig.add_subplot(grid[3, i], projection='polar')
            circle_ax.set_title('n=' + str(i) # + '\nE=' + str() + '$E_0$'
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
                        eigenstates[0:states].reshape([
                            states, N*N]))
            psi = psi.reshape([N, N])
            im.set_data(complex_to_rgb(psi))
            # apsi = np.abs(psi)
            # im.set_alpha(apsi/np.amax(apsi))
            # if animation_data['ticks'] % 2:
            #     return (im, )
            # else:
            for i, c in enumerate(coeffs):
                phi, r = np.angle(c), np.abs(c)
                artists[i].set_xdata([phi, phi])
                artists[i].set_ydata([0.0, r])
            return artists

        a = animation.FuncAnimation(fig, func, blit=True, interval=1000.0/60.0)
        plt.show()

