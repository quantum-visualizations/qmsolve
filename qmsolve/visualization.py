import numpy as np
import matplotlib.pyplot as plt


def visualize(energies, eigenstates, k):
    # 2d prototype for static visualization.  This should handle several the other dimensions

    fig = plt.figure(figsize=(16/9 *5.804 * 0.9,5.804)) 

    grid = plt.GridSpec(2, 2, width_ratios=[4.5, 1], height_ratios=[1, 1] , hspace=0.1, wspace=0.2)
    ax1 = fig.add_subplot(grid[0:2, 0:1])
    ax2 = fig.add_subplot(grid[0:2, 1:2])

    ax1.set_aspect('equal')
    ax1.set_ylabel("[Å]")
    ax1.set_xlabel("[Å]")

    ax2.set_title('E Level')
    ax2.set_facecolor('black')

    ax2.set_ylabel('$E_N$ (Relative to $E_{1}$)')
    ax2.set_xticks(ticks=[])


    E0 = energies[0]
    for E in energies:
        ax2.plot([0,1], [E/E0, E/E0], color='gray', alpha=0.5)

        
    im = ax1.imshow((eigenstates[k]), cmap ='seismic',  interpolation = 'bilinear')
    ax2.plot([0,1], [energies[k]/E0, energies[k]/E0], color='yellow', lw = 3)
    plt.show()





def dynamic_visualize(energies, eigenstates):
    # 2d slider prototype. This should handle several the other dimensions

    fig = plt.figure(figsize=(16/9 *5.804 * 0.9,5.804)) 

    grid = plt.GridSpec(2, 2, width_ratios=[5, 1], height_ratios=[1, 1] , hspace=0.1, wspace=0.2)
    ax1 = fig.add_subplot(grid[0:2, 0:1])
    ax2 = fig.add_subplot(grid[0:2, 1:2])

    ax1.set_aspect('equal')
    ax1.set_ylabel("[Å]")
    ax1.set_xlabel("[Å]")

    ax2.set_title('E Level')
    ax2.set_facecolor('black')

    ax2.set_ylabel('$E_N$ (Relative to $E_{1}$)')
    ax2.set_xticks(ticks=[])


    E0 = energies[0]
    for E in energies:
        ax2.plot([0,1], [E/E0, E/E0], color='gray', alpha=0.5)

    im = ax1.imshow((eigenstates[1]), cmap ='seismic',  interpolation = 'bilinear')
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
        im.set_data(eigenstates[state])
        line[0].set_ydata([energies[state]/E0, energies[state]/E0])

    slider.on_changed(update)
    plt.show()




def animate(energies, eigenstates):
    raise NotImplementedError('animate')