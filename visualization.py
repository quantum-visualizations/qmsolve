import numpy as np
import matplotlib.pyplot as plt

k = 3.8099821161548593 # hbar**2 / (2*m_e) /(Å**2) / eV
m_e = 1

def visualize(energies, eigenstates, k):
    # 2d prototype. This should be interactive and handle several the other dimensions

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
