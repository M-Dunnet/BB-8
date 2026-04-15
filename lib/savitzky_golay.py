#!/usr/bin/env python3

## Built-in importd
from math import factorial

## Installed imports
import numpy as np

def savitzky_golay(y, window_size, order, deriv=0, rate=1):
    '''
    Smooths over data using a Savitzky Golay filter
    This can either return a list of scores, or a list of peaks

    y:               Array-like, score list
    window_size:     INT, how big of a window to smooth
    order:           What order polynomial to use
    '''
    y = np.array(y)
    try:
        window_size = np.abs(int(window_size))
        order = np.abs(int(order))
    except ValueError:
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order + 1)
    half = (window_size - 1) // 2
    
    ## precompute coefficients
    b = np.array([[k**i for i in order_range] for k in range(-half, half + 1)])
    m = np.linalg.pinv(b)[deriv] * rate**deriv * factorial(deriv) 
    
    ## pad the signal at the extremes with values taken from the signal itself
    firstvals = y[0] - np.abs(y[1:half + 1][::-1] - y[0])
    lastvals = y[-1] + np.abs(y[-half - 1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    filtered = np.convolve(m[::-1], y, mode='valid')

    return filtered
