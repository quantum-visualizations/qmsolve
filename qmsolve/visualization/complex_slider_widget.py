import matplotlib.pyplot as plt
from matplotlib import widgets


class ComplexSliderWidget(widgets.AxesWidget):
    """
    A circular complex slider widget for manipulating complex
    values.

    References:
    - https://matplotlib.org/stable/api/widgets_api.
    - https://github.com/matplotlib/matplotlib/blob/
    1ba3ff1c273bf97a65e19892b23715d19c608ae5/lib/matplotlib/widgets.py
    """

    def __init__(self, ax, angle, r, animated=False):
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
            if not self.animated:
                event.canvas.draw()
