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


def complex_to_rgba(Z: np.ndarray, max_val: float = 1.0) -> np.ndarray:
    r = np.abs(Z)
    arg = np.angle(Z)
    
    h = (arg + np.pi)  / (2 * np.pi)
    s = np.ones(h.shape)
    v = np.ones(h.shape)  #alpha
    rgb = hsv_to_rgb(   np.moveaxis(np.array([h,s,v]) , 0, -1)  ) # --> tuple

    abs_z = np.abs(Z)/ max_val
    abs_z = np.where(abs_z> 1., 1. ,abs_z)
    return np.concatenate((rgb, abs_z.reshape((*abs_z.shape,1))), axis= (abs_z.ndim))
