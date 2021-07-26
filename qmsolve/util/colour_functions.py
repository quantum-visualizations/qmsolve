from matplotlib.colors import hsv_to_rgb
import numpy as np

def complex_to_rgb(Z):
    """Convert complex values to their rgb equivalent.

    Parameters
    ----------
    Z : array_like
        The complex values.

    Returns
    -------
    array_like
        The rgb values.
    """
    #using HSV space
    r = np.abs(Z)
    arg = np.angle(Z)
    
    h = (arg + np.pi)  / (2 * np.pi)
    s = np.ones(h.shape)
    v = r  / np.amax(r)  #alpha
    c = hsv_to_rgb(   np.moveaxis(np.array([h,s,v]) , 0, -1)  ) # --> tuple
    return c
