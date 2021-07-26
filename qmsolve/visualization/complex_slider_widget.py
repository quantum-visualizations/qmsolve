import matplotlib.pyplot as plt
from matplotlib import widgets


class ComplexSliderWidget(widgets.AxesWidget):
    """A circular complex slider widget for manipulating complex
    values.

    References:
    - https://matplotlib.org/stable/api/widgets_api.
    - https://github.com/matplotlib/matplotlib/blob/
    1ba3ff1c273bf97a65e19892b23715d19c608ae5/lib/matplotlib/widgets.py

    Parameters
    ----------
    widgets : matplotlib.widgets.AxesWidget
        A matplotlib widget.
    """

    def __init__(self, ax, angle, r, animated=False):
        """Constructor for the complex slider widgets.
        This takes the phase and magnitude of an initial complex number.

        Parameters
        ----------
        ax : matplotlib.axes
            plot object that holds the widget plot
        angle : float
            Phase of the initial complex number
        r : float
            Magnitude of the complex number
        animated : bool, optional
            Whether to animate this widget, by default False
        """
        line, = ax.plot([angle, angle], [0.0, r], linewidth=2.0)
        super().__init__(ax)
        self._rotator = line
        self._is_click = False
        self.animated = animated
        self.update = lambda x, y: None
        self.connect_event('button_press_event', self._click)
        self.connect_event('button_release_event', self._release)
        self.connect_event('motion_notify_event', self._motion)

    def get_artist(self):
        """Get the artist.

        Returns
        -------
        matplotlib.lines.Line2D
            Line plot whose length represents the magnitude
            of the complex value and its angular displacement
            represents its phase.
        """
        return self._rotator

    def _click(self, event):
        self._is_click = True
        self._update_plots(event)

    def _release(self, event):
        self._is_click = False

    def on_changed(self, update):
        """Set the update function which is used
        whenever the complex slider widget is changed.

        Parameters
        ----------
        update : callable
            The update function. It takes the phase
            and absolute magnitude of the complex value
            that this widget represents.
        """
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
            if not self.animated:
                event.canvas.draw()
