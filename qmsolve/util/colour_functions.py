from matplotlib.colors import hsv_to_rgb, to_rgba_array
import numpy as np

def complex_to_rgb(Z):
    #using HSV space
    r = np.abs(Z)
    arg = np.angle(Z)
    
    h = (arg + np.pi)  / (2 * np.pi)
    s = np.ones(h.shape)
    v = r  / np.amax(r)  #alpha
    c = hsv_to_rgb(   np.moveaxis(np.array([h,s,v]) , 0, -1)  ) # --> tuple
    return c


def complex_to_rgba(Z: np.ndarray, max_val: float = None) -> np.ndarray:
    alpha = ((np.abs(Z))/max_val if max_val else 
             np.abs(Z)/np.amax(np.abs(Z)))
    c1 = complex_to_rgb(Z)
    c2 = np.reshape(c1, (np.product(c1.shape[0:-1]), 3))
    a_flat = alpha.flatten()
    a_flat = np.where(a_flat <= 1.0, a_flat, np.exp(0.0*a_flat))
    c3 = to_rgba_array(c2, a_flat)
    new_shape = c1.shape[0:-1] + (4, )
    return c3.reshape(new_shape)
