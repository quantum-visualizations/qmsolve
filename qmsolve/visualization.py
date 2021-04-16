import numpy as np
import matplotlib.pyplot as plt


def visualize(energies, eigenstates, k):

    ndim = len(eigenstates[k].shape)
    # 2d prototype for static visualization.  This should handle several the other dimensions

    fig = plt.figure(figsize=(16/9 *5.804 * 0.9,5.804)) 

    grid = plt.GridSpec(2, 2, width_ratios=[4.5, 1], height_ratios=[1, 1] , hspace=0.1, wspace=0.2)
    ax1 = fig.add_subplot(grid[0:2, 0:1])
    ax2 = fig.add_subplot(grid[0:2, 1:2])

    ax1.set_xlabel("[Å]")

    ax2.set_title('E Level')
    ax2.set_facecolor('black')

    ax2.set_ylabel('$E_N$ (Relative to $E_{1}$)')
    ax2.set_xticks(ticks=[])


    E0 = energies[0]
    for E in energies:
        ax2.plot([0,1], [E/E0, E/E0], color='gray', alpha=0.5)

    if ndim == 2:
        ax1.set_aspect('equal')
        ax1.set_ylabel("[Å]")
        im = ax1.imshow((eigenstates[k]), cmap ='seismic',  interpolation = 'bilinear',  origin='lower')
    else:
        ax1.plot(eigenstates[k])

    ax2.plot([0,1], [energies[k]/E0, energies[k]/E0], color='yellow', lw = 3)
    plt.show()





def dynamic_visualize(energies, eigenstates):
    # 1d and 2d slider prototype. This should handle several the other dimensions

    ndim = len(eigenstates[0].shape)
    fig = plt.figure(figsize=(16/9 *5.804 * 0.9,5.804)) 

    grid = plt.GridSpec(2, 2, width_ratios=[5, 1], height_ratios=[1, 1] , hspace=0.1, wspace=0.2)
    ax1 = fig.add_subplot(grid[0:2, 0:1])
    ax2 = fig.add_subplot(grid[0:2, 1:2])

    ax1.set_xlabel("[Å]")

    ax2.set_title('E Level')
    ax2.set_facecolor('black')

    ax2.set_ylabel('$E_N$ (Relative to $E_{1}$)')
    ax2.set_xticks(ticks=[])


    E0 = energies[0]
    for E in energies:
        ax2.plot([0,1], [E/E0, E/E0], color='gray', alpha=0.5)
    
    if ndim == 2:
        ax1.set_aspect('equal')
        ax1.set_ylabel("[Å]")
        eigenstate_plot = ax1.imshow((eigenstates[1]), cmap ='seismic',  
                                     interpolation = 'bilinear', aspect ="equal",  extent = [-1,1,-1,1],  origin='lower' )
    elif ndim == 1:
        # ax1.set_xlim( ??? )
        eigenstate_plot, = ax1.plot(eigenstates[1])
        eigenstate_plot.set_data = eigenstate_plot.set_ydata

    line = ax2.plot([0,1], [energies[1]/E0, energies[1]/E0], color='yellow', lw = 3)

    plt.subplots_adjust(bottom=0.2)
    from matplotlib.widgets import Slider
    slider_ax = plt.axes([0.2, 0.05, 0.7, 0.05])
    slider = Slider(slider_ax,      # the axes object containing the slider
                      'state',            # the name of the slider parameter
                      0,          # minimal value of the parameter
                      len(eigenstates)-1,          # maximal value of the parameter
                      valinit = 1,  # initial value of the parameter 
                      valstep = 1,
                      color = '#5c05ff' 
                     )

    def update(state):
        state = int(state)
        eigenstate_plot.set_data(eigenstates[state])
        line[0].set_ydata([energies[state]/E0, energies[state]/E0])

    slider.on_changed(update)
    plt.show()




def animate(energies, eigenstates):

    ndim = len(eigenstates[0].shape)
    fig = plt.figure(figsize=(16/9 *5.804 * 0.9,5.804)) 

    grid = plt.GridSpec(2, 2, width_ratios=[5, 1], height_ratios=[1, 1] , hspace=0.1, wspace=0.2)
    ax1 = fig.add_subplot(grid[0:2, 0:1])
    ax2 = fig.add_subplot(grid[0:2, 1:2])

    ax1.set_xlabel("[Å]")

    ax2.set_title('E Level')
    ax2.set_facecolor('black')

    ax2.set_ylabel('$E_N$ (Relative to $E_{1}$)')
    ax2.set_xticks(ticks=[])


    E0 = energies[0]
    for E in energies:
        ax2.plot([0,1], [E/E0, E/E0], color='gray', alpha=0.5)
    
    if ndim == 2:
        ax1.set_aspect('equal')
        ax1.set_ylabel("[Å]")
        eigenstate_plot = ax1.imshow((eigenstates[1]), cmap ='seismic',  
                                     interpolation = 'bilinear',  origin='lower')
    elif ndim == 1:
        # ax1.set_xlim( ??? )
        eigenstate_plot, = ax1.plot(eigenstates[1])
        eigenstate_plot.set_data = eigenstate_plot.set_ydata

    line, = ax2.plot([0,1], [energies[1]/E0, energies[1]/E0], color='yellow', lw = 3)

    plt.subplots_adjust(bottom=0.2)

    import matplotlib.animation as animation

    animation_data = {'n': 0.0}
    def func_animation(*arg):
        animation_data['n'] = (animation_data['n'] + 0.1) % len(energies)
        state = int(animation_data['n'])
        if (animation_data['n'] % 1.0) > 0.5:
            transition_time = (animation_data['n'] - int(animation_data['n']) - 0.5)
            eigenstate_plot.set_data(np.cos(np.pi*transition_time)*eigenstates[state] + 
                                     np.sin(np.pi*transition_time)*
                                     eigenstates[(state + 1) % len(energies)])
            E_N = energies[state]/E0 
            E_M = energies[(state + 1) % len(energies)]/E0
            E =  E_N*np.cos(np.pi*transition_time)**2 + E_M*np.sin(np.pi*transition_time)**2
            line.set_ydata([E, E])
        else:
            line.set_ydata([energies[state]/E0, energies[state]/E0])
            eigenstate_plot.set_data(eigenstates[int(state)])
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


def compare_numerical_energies_with_analyitic(num_energies, analytic_energies):
    """
    Produce a plot that compares the numerical energies with 
    known analytical versions.

    Note that as of now this function is mostly a reuse from my plot 
    function from the earlier SPS-QM-2D project.
    """
    k = len(analytic_energies)
    fig = plt.figure()
    axes = fig.subplots(2, 1)
    axes[0].scatter([i for i in range(len(analytic_energies))],
                analytic_energies, marker='_', label='Analytic')
    axes[0].scatter([i for i in range(len(num_energies))], num_energies, marker='_')
    axes[0].set_xlim(-1, k)
    axes[0].set_ylim(0, np.amax(analytic_energies)*1.1)
    axes[0].set_ylabel('Energy')
    axes[0].grid()
    axes[0].legend()
    axes[1].set_xlim(-1, k)
    axes[1].scatter([i for i in range(len(num_energies))],
                     num_energies/np.array(analytic_energies), marker='+')
    axes[1].set_ylabel('numerical divided by analytic')
    axes[1].set_ylim(0.8, 1.2)
    axes[1].set_xlabel('Energy Count')
    axes[1].grid()
    plt.show()
    plt.close()