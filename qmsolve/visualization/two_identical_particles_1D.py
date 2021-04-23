import numpy as np
import matplotlib.pyplot as plt
from matplotlib import widgets
from matplotlib import animation
from .visualization import Visualization


from matplotlib.colors import hsv_to_rgb
def complex_to_rgb(Z):
    #using HSV space
    r = np.abs(Z)
    arg = np.angle(Z)
    
    h = (arg + np.pi)  / (2 * np.pi)
    s = np.ones(h.shape)
    v = r  / np.amax(r)  #alpha
    c = hsv_to_rgb(   np.moveaxis(np.array([h,s,v]) , 0, -1)  ) # --> tuple
    return c



class VisualizationIdenticalParticles1D(Visualization):
    def __init__(self,eigenstates):
        self.eigenstates = eigenstates



    def plot_eigenstate(self, k):
        eigenstates_array = self.eigenstates.array
        energies = self.eigenstates.energies
        plt.style.use("dark_background")

        fig = plt.figure(figsize=(7.5 ,7.0)) 

        grid = plt.GridSpec(2, 2, width_ratios=[4.5, 1], height_ratios=[2.5, 1] , hspace=0.4, wspace=0.4)
        ax1 = fig.add_subplot(grid[0:1, 0:1])
        ax3 = fig.add_subplot(grid[1:2, 0:1], sharex=ax1) # probability density of finding any particle at x 
        ax2 = fig.add_subplot(grid[0:2, 1:2])

        ax1.set_xlabel("$x_1$ [Å]")
        ax1.set_ylabel("$x_2$ [Å]")
        ax1.set_title("$\Psi(x_1,x_2)$")

        ax2.set_title('E Level')
        ax2.set_facecolor('black')
        ax2.set_ylabel('$E_N$ (Relative to $E_{1}$)')
        ax2.set_xticks(ticks=[])

        ax3.set_xlabel("$x$ [Å]")
        ax3.set_ylabel("${\| \Psi(x)\|}^{2} $")
        ax3.set_title("Probability density")

        plt.setp(ax3.get_yticklabels(), visible=False)


        E0 = energies[0]

        for E in energies:
            ax2.plot([0,1], [E/E0, E/E0], color='gray', alpha=0.5)

        ax2.plot([0,1], [energies[k]/E0, energies[k]/E0], color='yellow', lw = 3)


        L =  self.eigenstates.extent/2
        im = ax1.imshow(complex_to_rgb(eigenstates_array[k]*np.exp( 1j*2*np.pi/10*k)), aspect = "auto", origin='lower',extent = [-L, L, -L, L],  interpolation = 'bilinear')
        x = np.linspace(-L, L, self.eigenstates.N)

        prob_density = np.abs(np.sum(  (eigenstates_array[k])*np.conjugate(eigenstates_array[k])  , axis = 1))
        prob_plot = ax3.plot(x,  prob_density, color= "cyan")
        prob_plot_fill = ax3.fill_between(x,prob_density, alpha=0.1, color= "cyan" )
        ax3.set_ylim([0,max(prob_density)*1.1])
        plt.show()


    def slider_plot(self):

        eigenstates_array = self.eigenstates.array
        energies = self.eigenstates.energies
        plt.style.use("dark_background")

        fig = plt.figure(figsize=(7.5 ,7.0)) 

        grid = plt.GridSpec(2, 2, width_ratios=[4.5, 1], height_ratios=[2.5, 1] , hspace=0.4, wspace=0.4)
        ax1 = fig.add_subplot(grid[0:1, 0:1])
        ax3 = fig.add_subplot(grid[1:2, 0:1], sharex=ax1) # probability density of finding any particle at x 
        ax2 = fig.add_subplot(grid[0:2, 1:2])

        ax1.set_xlabel("$x_1$ [Å]")
        ax1.set_ylabel("$x_2$ [Å]")
        ax1.set_title("$\Psi(x_1,x_2)$")

        ax2.set_title('E Level')
        ax2.set_facecolor('black')
        ax2.set_ylabel('$E_N$ (Relative to $E_{1}$)')
        ax2.set_xticks(ticks=[])

        ax3.set_xlabel("$x$ [Å]")
        ax3.set_ylabel("${\| \Psi(x)\|}^{2} $")
        ax3.set_title("Probability density")
        plt.setp(ax3.get_yticklabels(), visible=False)

        E0 = energies[0]
        for E in energies:
            ax2.plot([0,1], [E/E0, E/E0], color='gray', alpha=0.5)


        L =  self.eigenstates.extent/2
        eigenstate_plot = ax1.imshow(complex_to_rgb(eigenstates_array[0]*np.exp( 1j*2*np.pi/10*0)), aspect = "auto", origin='lower',extent = [-L, L, -L, L],   interpolation = 'bilinear')
        x = np.linspace(-L, L, self.eigenstates.N)

        prob_density = np.abs(np.sum(  (eigenstates_array[0])*np.conjugate(eigenstates_array[0])  , axis = 1))
        prob_plot = ax3.plot(x,  prob_density, color= "cyan")
        prob_plot_fill = ax3.fill_between(x,prob_density, alpha=0.1, color= "cyan" )
        
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


            ax3.clear()
            ax3.set_xlabel("$x$ [Å]")
            ax3.set_ylabel("${\| \Psi(x)\|}^{2} $")
            ax3.set_title("Probability density")
            ax3.set_yticks([])
            plt.setp(ax3.get_yticklabels(), visible=False)

            prob_density = np.abs(np.sum(  (eigenstates_array[state])*np.conjugate(eigenstates_array[state])  , axis = 1))
            prob_plot = ax3.plot(x,  prob_density, color= "cyan")
            ax3.set_ylim([0,max(prob_density)*1.1])
            prob_plot_fill = ax3.fill_between(x,prob_density, alpha=0.1, color= "cyan" )
            line[0].set_ydata([energies[state]/E0, energies[state]/E0])
        slider.on_changed(update)
        plt.show()






    def animate(self, max_states = None):

        print(len(self.eigenstates.energies))
        if max_states == None:
            max_states = len(self.eigenstates.energies)

        eigenstates_array = self.eigenstates.array[:max_states]
        energies = self.eigenstates.energies[:max_states]

        plt.style.use("dark_background")

        fig = plt.figure(figsize=(7.5 ,7.0)) 

        grid = plt.GridSpec(2, 2, width_ratios=[4.5, 1], height_ratios=[2.5, 1] , hspace=0.4, wspace=0.4)
        ax1 = fig.add_subplot(grid[0:1, 0:1])
        ax3 = fig.add_subplot(grid[1:2, 0:1], sharex=ax1) # probability density of finding any particle at x 
        ax2 = fig.add_subplot(grid[0:2, 1:2])

        ax1.set_xlabel("$x_1$ [Å]")
        ax1.set_ylabel("$x_2$ [Å]")
        ax1.set_title("$\Psi(x_1,x_2)$")

        ax2.set_title('E Level')
        ax2.set_facecolor('black')
        ax2.set_ylabel('$E_N$ (Relative to $E_{1}$)')
        ax2.set_xticks(ticks=[])

        ax3.set_xlabel("$x$ [Å]")
        ax3.set_ylabel("${\| \Psi(x)\|}^{2} $")
        ax3.set_title("Probability density")
        #ax3.set_yticks([])
        plt.setp(ax3.get_yticklabels(), visible=False)

        E0 = energies[0]
        for E in energies:
            ax2.plot([0,1], [E/E0, E/E0], color='gray', alpha=0.5)


        L =  self.eigenstates.extent/2
        eigenstate_plot = ax1.imshow(complex_to_rgb(eigenstates_array[0]), aspect = "auto",origin='lower', extent = [-L, L, -L, L],   interpolation = 'bilinear')
        x = np.linspace(-L, L, self.eigenstates.N)

        #prob_density = np.sum(  (eigenstates_array[0])*np.conjugate(eigenstates_array[0])  , axis = 1)
        #ax3.set_ylim([0,max(prob_density)*1.2])
        #prob_plot, = ax3.plot(x,  prob_density, color= "cyan")
        #prob_plot_fill = ax3.fill_between(x,prob_density, alpha=0.1, color= "cyan" )
        line, = ax2.plot([0,1], [energies[0]/E0, energies[0]/E0], color='yellow', lw = 3)

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

                prob_density = np.abs(np.sum(  eigenstate_combination*np.conjugate(eigenstate_combination)  , axis = 1))
                prob_plot, = ax3.plot(x,  prob_density, color= "cyan")
                prob_plot_fill = ax3.fill_between(x,prob_density, alpha=0.1, color= "cyan" )
                ax3.set_ylim([0,max(prob_density)*1.1])
                E_N = energies[state]/E0 
                E_M = energies[(state + 1) % len(energies)]/E0
                E =  E_N*np.cos(np.pi*transition_time)**2 + E_M*np.sin(np.pi*transition_time)**2
                line.set_ydata([E, E])
            else:
                line.set_ydata([energies[state]/E0, energies[state]/E0])

                eigenstate_combination = eigenstates_array[int(state)]*np.exp( 1j*2*np.pi/10*state)
                eigenstate_plot.set_data(complex_to_rgb(eigenstate_combination))
                
                prob_density = np.abs(np.sum(  eigenstate_combination*np.conjugate(eigenstate_combination)  , axis = 1))
                prob_plot, = ax3.plot(x,  prob_density, color= "cyan")
                prob_plot_fill = ax3.fill_between(x,prob_density, alpha=0.1, color= "cyan" )
                ax3.set_ylim([0,max(prob_density)*1.1])
            return eigenstate_plot, line,prob_plot,prob_plot_fill

        a = animation.FuncAnimation(fig, func_animation,
                                    blit=True, interval=1.0)
        plt.show()
        
        # save animation
        """
        Writer = animation.writers['ffmpeg']
        writer = Writer(fps=20, metadata=dict(artist='Me'), bitrate=1800)
        a.save('im.mp4', writer=writer)
        """        
