#!/usr/bin/env python3
'''
Logic for calling repeated splint peaks. Peaks determine where the long ONT reads will
be split into subreads.
'''
## BB-8 imports
from lib.savitzky_golay import savitzky_golay

## Installed imports
from scipy.signal import find_peaks
import numpy as np

def call_peaks(scores, min_dist, iters, window, order):
    '''
    Identify peaks in a 1D score array after iterative Savitzky Golay smoothing.
    
    Scores are define by conk. The peak-finder argument sets the parameters for this function:
    Deafults to: 20,3,41,2 | min_dist=20, iters=3, window=41, order=2

    The signal is smoothed `iters` times, then peaks are detected using a
    minimum distance (`min_dist`) and a height threshold based on the median
    score. Returns no peaks if the signal lacks sufficient prominence.
    '''
    peaks = []
    for i in range(iters):
        scores = savitzky_golay(scores, window, order, deriv=0, rate=1)
    med_score = np.median(scores)
    if max(scores) < 6 * med_score:
        return peaks
    peaks, _ = find_peaks(scores, distance=min_dist, height=med_score * 3)
    return peaks
