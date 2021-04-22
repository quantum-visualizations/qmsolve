import numpy as np
import matplotlib.pyplot as plt
from matplotlib import widgets
from matplotlib import animation


class CircleWidget(widgets.AxesWidget):
    """
    Widget for modifying a complex value.
    """

    def __init__(self, ax, angle, r):
        line, = ax.plot([angle, angle], [0.0, r], linewidth=2.0)
        super().__init__(ax)
        self._rotator = line
        self._is_click = False
        self.update = lambda x, y: None
        self.connect_event('button_press_event', self._click)
        self.connect_event('button_release_event', self._release)
        self.connect_event('motion_notify_event', self._motion)

    def get_artist(self):
        return self._rotator

    def _click(self, event):
        self._is_click = True
        self._update_plots(event)

    def _release(self, event):
        self._is_click = False

    def on_changed(self, update):
        self.update = update
    
    def _motion(self, event):
        self._update_plots(event)

    def _update_plots(self, event):
        if (self._is_click and event.xdata != None
            and event.ydata != None
            and event.x >= self.ax.bbox.xmin and
            event.x < self.ax.bbox.xmax and
            event.y >= self.ax.bbox.ymin and
            event.y < self.ax.bbox.ymax
            ):
            phi, r = event.xdata, event.ydata 
            if r < 0.2:
                r = 0.0
            self.update(phi, r)
            self._rotator.set_xdata([phi, phi])
            self._rotator.set_ydata([0.0, r])


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



def visualize_superpositions(energies, eigenstates, n_states, **kw):
    """
    Visualize the time evolution of a superposition of energy eigenstates.
    The circle widgets control the relative phase of each of the eigenstates.
    These widgets are inspired by the circular phasors from the
    quantum mechanics applets by Paul Falstad:
    https://www.falstad.com/qm1d/

    """
    params = {'dt': 0.001}
    for k in kw.keys():
        params[k] = kw[k]
    ndim = len(eigenstates[0].shape)
    if ndim == 2:
        _visualize_superpositions_2d(energies, eigenstates, n_states, params)
        return
    fig = plt.figure(figsize=(16/9 *5.804 * 0.9,5.804)) 
    grid = plt.GridSpec(4, 10)
    ax = fig.add_subplot(grid[0:3, 0:10])
    ax.set_xticks([])
    ax.set_yticks([])
    get_norm_factor = lambda psi: 1.0/np.sqrt(np.sum(psi*np.conj(psi)))
    coeffs = np.array([1.0 if i == 0 else 0.0 for i in range(n_states)],
                      dtype=np.complex128)
    ax.set_ylim(-1.7*np.amax(eigenstates[0]), 
                1.7*np.amax(eigenstates[0]))
    line1, = ax.plot(np.real(eigenstates[0]))
    line2, = ax.plot(np.imag(eigenstates[0]))
    line3, = ax.plot(np.abs(eigenstates[0]), color='black')
    animation_data = {'ticks': 0, 'norm': 1.0}

    def make_update(n):
        def update(phi, r):
            coeffs[n] = r*np.exp(1.0j*phi)
            psi = np.tensordot(coeffs, eigenstates[0:n_states], 1)
            animation_data['norm'] = get_norm_factor(psi)
            psi *= animation_data['norm']
            line1.set_ydata(np.real(psi))
            line2.set_ydata(np.imag(psi))
            line3.set_ydata(np.abs(psi))
        return update

    widgets = []
    circle_artists = []
    for i in range(n_states):
        circle_ax = fig.add_subplot(grid[3, i], projection='polar')
        circle_ax.set_title('n=' + str(i) # + '\nE=' + str() + '$E_0$'
                            )
        circle_ax.set_xticks([])
        circle_ax.set_yticks([])
        widgets.append(CircleWidget(circle_ax, 0.0, 1.0))
        widgets[i].on_changed(make_update(i))
        circle_artists.append(widgets[i].get_artist())
    artists = circle_artists + [line1, line2, line3]

    def func(*args):
        animation_data['ticks'] += 1
        e = np.exp(-1.0j*energies[0:n_states]*params['dt'])
        np.copyto(coeffs, coeffs*e)
        norm_factor = animation_data['norm']
        psi = 4*np.tensordot(coeffs*norm_factor, eigenstates[0:n_states], 1)
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

    a = animation.FuncAnimation(fig, func, blit=True, interval=1000.0/60.0)
    plt.show()



def _visualize_superpositions_2d(energies, eigenstates, n_states, params):
    eigenstates = np.array(eigenstates)
    energies = np.array(energies)
    N = eigenstates.shape[1]
    fig = plt.figure(figsize=(16/9 *5.804 * 0.9,5.804)) 
    grid = plt.GridSpec(4, 10)
    ax = fig.add_subplot(grid[0:3, 0:10])
    ax.set_xticks([])
    ax.set_yticks([])
    get_norm_factor = lambda psi: 1.0/np.sqrt(np.sum(psi*np.conj(psi)))
    coeffs = np.array([1.0 if i == 0 else 0.0 for i in range(n_states)],
                    dtype=np.complex128)
    X, Y = np.meshgrid(np.linspace(-1.0, 1.0, eigenstates[0].shape[0]),
                       np.linspace(-1.0, 1.0, eigenstates[0].shape[1]))
    maxval = np.amax(np.abs(eigenstates[0]))
    im = plt.imshow(np.angle(X + 1.0j*Y), alpha=np.abs(eigenstates[0])/maxval, 
                    cmap='hsv', interpolation='none')
    im2 = plt.imshow(0.0*eigenstates[0], cmap='gray')
    animation_data = {'ticks': 0, 'norm': 1.0}

    def make_update(n):
        def update(phi, r):
            coeffs[n] = r*np.exp(1.0j*phi)
            psi = np.dot(coeffs, 
                         eigenstates[0:n_states].reshape([n_states, N*N]))
            psi = psi.reshape([N, N])
            animation_data['norm'] = get_norm_factor(psi)
            psi *= animation_data['norm']
            # apsi = np.abs(psi)
            # im.set_alpha(apsi/np.amax(apsi))
        return update

    widgets = []
    circle_artists = []
    for i in range(n_states):
        circle_ax = fig.add_subplot(grid[3, i], projection='polar')
        circle_ax.set_title('n=' + str(i) # + '\nE=' + str() + '$E_0$'
                            )
        circle_ax.set_xticks([])
        circle_ax.set_yticks([])
        widgets.append(CircleWidget(circle_ax, 0.0, 1.0))
        widgets[i].on_changed(make_update(i))
        circle_artists.append(widgets[i].get_artist())
    artists = circle_artists + [im]

    def func(*args):
        animation_data['ticks'] += 1
        e = np.exp(-1.0j*energies[0:n_states]*params['dt'])
        np.copyto(coeffs, coeffs*e)
        norm_factor = animation_data['norm']
        psi = np.dot(coeffs*norm_factor, 
                     eigenstates[0:n_states].reshape([n_states, N*N]))
        psi = psi.reshape([N, N])
        im.set_data(np.angle(psi))
        apsi = np.abs(psi)
        im.set_alpha(apsi/np.amax(apsi))
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


def compare_numerical_with_analyitic(num_energies, analytic_energies):
    """
    Produce a plot that compares the numerical energies with 
    known analytical versions.
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